package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.core.java.io.LoadNanoJTable;
import nanoj.core2.NanoJProfiler;

import java.io.IOException;
import java.nio.FloatBuffer;
import java.nio.IntBuffer;
import java.util.Map;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.min;
import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.getResourceAsString;
import static nanoj.core2.NanoJCL.replaceFirst;

public class liveSRRF_CL {

    // Basic formats
    private int width,
            height,
            widthM,
            heightM,
            nFrameForSRRF,
            blockLength,
            magnification,
            nPlanesThreeD,
            nSplits = 3, // TODO: currently hard-coded
            widthS,
            heightS;

    private double[] intCoeffsThreeD,
            chosenROIsLocationsThreeD;

    private final int GxGyMagnification = 2;
    private final float vxy_offset = 0.5f;
    private final int vxy_ArrayShift = 1;

    private final int nReconstructions = 2; // Currently only STD and AVG

    private boolean doMPmapCorrection,
            do3DSRRF;

    public ImageStack imsSRRF;

    // Advanced formats
    private final NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programLiveSRRF;
    static private CLKernel kernelCalculateGradient,
            kernelInterpolateGradient,
            kernelIncrementFramePosition,
            kernelResetFramePosition,
            kernelCalculateSRRF,
            kernelCalculateMPmap,
            kernelCalculateStd,
            kernelCorrectMPmap;

    static private CLPlatform clPlatformMaxFlop;
    static private CLDevice clDeviceMaxFlop;
    static public CLDevice[] allCLdevices;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy, clBufferGz,
            clBufferGxInt, clBufferGyInt, clBufferGzInt,
            clBufferShiftXY,
            clBufferShiftXYthreeD,
            clBufferOut,
            clBufferMPmap;

    private CLBuffer<IntBuffer>
            clBufferCurrentFrame;


    // --- Constructor ---
    public liveSRRF_CL() {
        // Nothing to see here. Keep calm and carry on.
    }


    // -- Check devices --
    public void checkDevices() {

        int nDevices = 0;
        CLPlatform[] allPlatforms;

        try{allPlatforms = CLPlatform.listCLPlatforms();}
        catch(CLException ex) {
            IJ.log("Something went wrong initializing OpenCL.");
            throw new RuntimeException("Something went wrong initializing OpenCL.");
        }

        double nFlops = 0;

        for (CLPlatform allPlatform : allPlatforms) {
            CLDevice[] allCLdeviceOnThisPlatform = allPlatform.listCLDevices();
            nDevices += allCLdeviceOnThisPlatform.length;

            for (CLDevice clDevice : allCLdeviceOnThisPlatform) {
                IJ.log("--------");
                IJ.log("Device name: " + clDevice.getName());
                IJ.log("Device type: " + clDevice.getType());
                IJ.log("Max clock: " + clDevice.getMaxClockFrequency() + " MHz");
                IJ.log("Number of compute units: " + clDevice.getMaxComputeUnits());
                if (clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency() > nFlops) {
                    nFlops = clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency();
                    clPlatformMaxFlop = allPlatform;
                    clDeviceMaxFlop = clDevice;
                }
            }
        }

        IJ.log("--------");
        IJ.log("Maximum flops device: " + clDeviceMaxFlop.getName());

        allCLdevices = new CLDevice[nDevices];

        int i = 0;
        for (CLPlatform allPlatform : allPlatforms) {
            CLDevice[] allCLdeviceOnThisPlatform = allPlatform.listCLDevices();
            nDevices += allCLdeviceOnThisPlatform.length;

            for (CLDevice clDevice : allCLdeviceOnThisPlatform) {
                allCLdevices[i] = clDevice;
                i++;
            }
        }

    }


    // --- Initialization method ---
    public void initialise(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF, int blockLength, CLDevice chosenDevice, boolean intWeighting, boolean doMPmapCorrection, String threeDimSRRFcalibTablePath) {

        this.width = width;
        this.height = height;
        this.heightM = height * magnification;
        this.widthM = width * magnification;
        this.nFrameForSRRF = nFrameForSRRF;
        this.blockLength = blockLength;
        this.magnification = magnification;
        this.doMPmapCorrection = doMPmapCorrection;

        if (chosenDevice == null) {
//            IJ.log("Looking for the fastest device...");
            System.out.println("Using the fastest device...");
            context = CLContext.create(clPlatformMaxFlop);
//            context = CLContext.create();
            chosenDevice = context.getMaxFlopsDevice();
//            IJ.log("Using "+chosenDevice.getName());
        }
        else{
            context = CLContext.create(chosenDevice.getPlatform());
            IJ.log("Working on platform: "+chosenDevice.getPlatform().getName());
            CLDevice[] allCLdevicesOnThisPlatform = context.getDevices();
            int i = 0;
            while (!allCLdevicesOnThisPlatform[i].getName().equals(chosenDevice.getName())){
                i++;
            }
            chosenDevice = allCLdevicesOnThisPlatform[i];
        }

        // 3D-SRRF initialisation
        do3DSRRF = threeDimSRRFcalibTablePath != null;

        if(do3DSRRF){
            double[] shiftXthreeD = null;
            double[] shiftYthreeD;
            double[] thetaThreeD;
            double[] axialPositionsThreeD;

            Map<String, double[]> calibTable;
            try {
                calibTable = new LoadNanoJTable(threeDimSRRFcalibTablePath).getData();
                shiftXthreeD = calibTable.get("X-shift (pixels)");
                shiftYthreeD = calibTable.get("Y-shift (pixels)");
                thetaThreeD = calibTable.get("Theta (degrees)");
                chosenROIsLocationsThreeD = calibTable.get("ROI #");
                axialPositionsThreeD = calibTable.get("Axial positions");
                intCoeffsThreeD = calibTable.get("Intensity scaling");
                ResultsTable rt = dataMapToResultsTable(calibTable);
                rt.show("Calibration-Table");
            } catch (IOException e) {
                IJ.log("Catching exception...");
                e.printStackTrace();
            }
            nPlanesThreeD = shiftXthreeD.length;
            widthS = width/nSplits;
            heightS = height/nSplits;
        }


        System.out.println("using " + chosenDevice);
        //IJ.log("Using " + chosenDevice.getName());

        // initialise buffers
        clBufferPx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_ONLY);
        clBufferShiftXY = context.createFloatBuffer(2 * nFrameForSRRF, READ_ONLY);
        clBufferGx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // all frames Gx
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // all frames Gy
        clBufferGxInt = context.createFloatBuffer(nFramesOnGPU * GxGyMagnification * GxGyMagnification * width * height, READ_WRITE); // single frame Gx
        clBufferGyInt = context.createFloatBuffer(nFramesOnGPU * GxGyMagnification * GxGyMagnification * width * height, READ_WRITE); // single frame Gy

        // initialise buffers for 3D-SRRF
        if(do3DSRRF) {
            clBufferShiftXYthreeD = context.createFloatBuffer(2 * nPlanesThreeD * nFrameForSRRF, READ_ONLY);
            // TODO: add some buffers for the rotation
            clBufferGz = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE);
            clBufferGzInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE);
            // TODO: this assumes only interpolation in XY, not in Z yet
        }

        clBufferOut = context.createFloatBuffer((nReconstructions + 1) * widthM * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);
        clBufferMPmap = context.createFloatBuffer(2 * magnification * magnification, READ_WRITE);

        // Current frame is a 2 element Int buffer:
        // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
        // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every tun of the method calculateSRRF (within the gradient calculation))

        // Create the program
        float sigma = fwhm / 2.354f;
        float radius = ((float) ((int) (GxGyMagnification * 2 * sigma))) / GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I

        String programString;
        if(do3DSRRF) programString = getResourceAsString(liveSRRF_CL.class, "live3DSRRF.cl");
        else programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");

        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
        programString = replaceFirst(programString, "$FWHM$", "" + fwhm);
        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
        programString = replaceFirst(programString, "$GXGYMAGNIFICATION$", "" + GxGyMagnification);

        programString = replaceFirst(programString, "$SIGMA$", "" + sigma);
        programString = replaceFirst(programString, "$RADIUS$", "" + radius);
        programString = replaceFirst(programString, "$WIDTH$", "" + width);
        programString = replaceFirst(programString, "$HEIGHT$", "" + height);
        programString = replaceFirst(programString, "$WH$", "" + (width * height));
        programString = replaceFirst(programString, "$WM$", "" + (width * magnification));
        programString = replaceFirst(programString, "$HM$", "" + (height * magnification));
        programString = replaceFirst(programString, "$WHM$", "" + (width * height * magnification * magnification));
        programString = replaceFirst(programString, "$WINT$", "" + (GxGyMagnification * width));
        programString = replaceFirst(programString, "$HINT$", "" + (GxGyMagnification * height));
        programString = replaceFirst(programString, "$VXY_OFFSET$", "" + vxy_offset);
        programString = replaceFirst(programString, "$VXY_ARRAYSHIFT$", "" + vxy_ArrayShift);
        programString = replaceFirst(programString, "$NFRAMEFORSRRF$", "" + nFrameForSRRF);

        if(do3DSRRF) {
            programString = replaceFirst(programString, "$WIDTHTHREED$", "" + width/nSplits); // TODO: split is hardcoded at the mo, 3x3 for MFM
            programString = replaceFirst(programString, "$HEIGHTTHREED$", "" + height/nSplits);
            programString = replaceFirst(programString, "$WHS$", "" + ((width/nSplits) * (height/nSplits)));
            programString = replaceFirst(programString, "$WSINT$", "" + ((width/nSplits)*GxGyMagnification));
            programString = replaceFirst(programString, "$HSINT$", "" + ((height/nSplits)*GxGyMagnification));
            programString = replaceFirst(programString, "$WHSINT$", "" + ((width/nSplits)*GxGyMagnification * (height/nSplits)*GxGyMagnification));
            programString = replaceFirst(programString, "$NPLANES$", "" + nPlanesThreeD);
        }

        if (intWeighting) programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 1);
        else programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 0);

        programLiveSRRF = context.createProgram(programString).build();

        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");
        kernelInterpolateGradient = programLiveSRRF.createCLKernel("calculateGradientInterpolation");
        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateRadialGradientConvergence");
        kernelIncrementFramePosition = programLiveSRRF.createCLKernel("kernelIncrementFramePosition");
        kernelResetFramePosition = programLiveSRRF.createCLKernel("kernelResetFramePosition");
        kernelCalculateMPmap = programLiveSRRF.createCLKernel("kernelCalculateMPmap");
        kernelCalculateStd = programLiveSRRF.createCLKernel("kernelCalculateStd");
        kernelCorrectMPmap = programLiveSRRF.createCLKernel("kernelCorrectMPmap");

        int argn;
        argn = 0;
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        if(do3DSRRF) kernelCalculateGradient.setArg(argn++, clBufferGz); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelInterpolateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        if(do3DSRRF) kernelInterpolateGradient.setArg(argn++, clBufferGz); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        if(do3DSRRF) kernelInterpolateGradient.setArg(argn++, clBufferGzInt); // make sure type is the same !!

        argn = 0;
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        if(do3DSRRF) kernelCalculateSRRF.setArg(argn++, clBufferGzInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftXY); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelIncrementFramePosition.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelCalculateMPmap.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateMPmap.setArg(argn++, clBufferMPmap); // make sure type is the same !!

        argn = 0;
        kernelCalculateStd.setArg(argn++, clBufferOut);

        argn = 0;
        kernelCorrectMPmap.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCorrectMPmap.setArg(argn++, clBufferMPmap); // make sure type is the same !!


        queue = chosenDevice.createCommandQueue();

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferShiftXY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferGxInt.getCLSize() +
                        clBufferGyInt.getCLSize() +
                        clBufferOut.getCLSize() +
                        clBufferCurrentFrame.getCLSize() +
                        clBufferMPmap.getCLSize()
        )
                / 1000000d + "MB");

    }


    // --- Load Shift array on GPU ---
    public void loadShiftXYGPUbuffer(float[] shiftX, float[] shiftY) {

        float[] shiftXY = new float[2 * nFrameForSRRF];
        System.arraycopy(shiftX, 0, shiftXY, 0, nFrameForSRRF);
        System.arraycopy(shiftY, 0, shiftXY, nFrameForSRRF, nFrameForSRRF);

        int id = prof.startTimer();
        fillBuffer(clBufferShiftXY, shiftXY);
        queue.putWriteBuffer(clBufferShiftXY, false);
        prof.recordTime("Uploading shift arrays to GPU", prof.endTimer(id));
    }


    // --- Calculate SRRF images ---
    public synchronized boolean calculateSRRF(ImageStack imsRawData) { // returns boolean describing whether it was cancelled by user or not

        assert (imsRawData.getWidth() == width && imsRawData.getHeight() == height);
        int nFrameToLoad = imsRawData.getSize();

        //        IJ.log("Uploading raw data to GPU...");
        int id = prof.startTimer();
        if(do3DSRRF){
            ImageStack imsData3D = reshapeImageStack3D(imsRawData);
            fillBuffer(clBufferPx, imsData3D);
        }
        else fillBuffer(clBufferPx, imsRawData);

        queue.putWriteBuffer(clBufferPx, false);
        prof.recordTime("Uploading data to GPU", prof.endTimer(id));

        // Make kernelCalculateGradient assignment (this kernel also resets the local GPU load frame counter)
//        IJ.log("Calculating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        if(do3DSRRF) queue.put1DRangeKernel(kernelCalculateGradient, 0,  width * height * nFrameToLoad, 0);
        else         queue.put3DRangeKernel(kernelCalculateGradient, 0, 0, 0, width, height, nFrameToLoad, 0, 0, 0);

        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

//        IJ.log("Interpolating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        if(do3DSRRF) queue.put1DRangeKernel(kernelInterpolateGradient, 0,  GxGyMagnification * width * GxGyMagnification * height * nFrameToLoad, 0);
        else         queue.put3DRangeKernel(kernelInterpolateGradient, 0, 0, 0, GxGyMagnification * width, GxGyMagnification * height, nFrameToLoad, 0, 0, 0);

        prof.recordTime("kernelInterpolateGradient", prof.endTimer(id));

        // Make kernelCalculateSRRF assignment
//        IJ.log("Calculating SRRF...");

        int workSize;
        int nBlocks = widthM * heightM / blockLength + ((widthM * heightM % blockLength == 0) ? 0 : 1);

        for (int f = 0; f < nFrameToLoad; f++) {

            for (int nB = 0; nB < nBlocks; nB++) {
                workSize = min(blockLength, widthM * heightM - nB * blockLength);

                id = prof.startTimer();
                queue.put1DRangeKernel(kernelCalculateSRRF, nB * blockLength, workSize, 0);
                prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return true;
                }

            }

            // This kernel needs to be done outside of the previous kernel because of concommitent execution (you never know when each pixel is executed)
            id = prof.startTimer();
            queue.finish(); // Make sure everything is done
            queue.put1DRangeKernel(kernelIncrementFramePosition, 0, 2, 0);
            prof.recordTime("Increment frame count", prof.endTimer(id));
        }

        return false;

    }

    // --- Read the output buffer ---
    public void readSRRFbuffer() {

        // Calculate the STD on the OutputArray on the GPU
        int id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put1DRangeKernel(kernelCalculateStd, 0, heightM*widthM,0);
        prof.recordTime("Calculate STD image", prof.endTimer(id));

        if (doMPmapCorrection) {
            id = prof.startTimer();
            queue.put1DRangeKernel(kernelCalculateMPmap, 0, 2 * magnification * magnification, 0);
            prof.recordTime("Calculate MP map", prof.endTimer(id));

            id = prof.startTimer();
            queue.put1DRangeKernel(kernelCorrectMPmap, 0, 2 * heightM*widthM, 0);
            prof.recordTime("Correct for MP map", prof.endTimer(id));
        }

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferOut, true);
        FloatBuffer bufferSRRF = clBufferOut.getBuffer();

        imsSRRF = new ImageStack(widthM, heightM);

        // Load average
        float[] dataSRRF = new float[widthM * heightM];
        for (int n = 0; n < widthM * heightM; n++) {
            dataSRRF[n] = bufferSRRF.get(n);
            if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we dont get any weirdness
        }
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));

        // Load standard deviation // TODO: do this on GPU - test but this looks fine
        dataSRRF = new float[widthM * heightM];
        for (int n = 0; n < widthM * heightM; n++) {
            //dataSRRF[n] = bufferSRRF.get(n + widthM * heightM) - bufferSRRF.get(n) * bufferSRRF.get(n); // Var[X] = E[X^2] - (E[X])^2
            dataSRRF[n] = bufferSRRF.get(n + widthM * heightM);
            //if (dataSRRF[n] < 0) {
            //    dataSRRF[n] = 0;
            //    //IJ.log("!!Negative VAR value!!");
            //}
            if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we dont get any weirdness
            //dataSRRF[n] = (float) math.sqrt(dataSRRF[n]);
        }
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));

        // Load interpolated data
        dataSRRF = new float[widthM * heightM];
        for (int n = 0; n < widthM * heightM; n++) {
            dataSRRF[n] = bufferSRRF.get(n + 2 * widthM * heightM);
            if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we don't get any weirdness
        }
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));

    }

    // --- Read the gradient buffers --- only used for testing!
    public ImageStack readGradientBuffers(boolean interpolated) {
        // order defines whether the output is the gradient or the interpolated gradients

        queue.finish(); // Make sure everything is done

        int imageWidth;
        int imageHeight;
        if(do3DSRRF){
            imageWidth = widthS;
            imageHeight = heightS;
        }
        else{
            imageWidth = width;
            imageHeight = height;
        }

        FloatBuffer bufferGx;
        FloatBuffer bufferGy;
        FloatBuffer bufferGz = null;

        if (!interpolated){
            queue.putReadBuffer(clBufferGx, true);
            queue.putReadBuffer(clBufferGy, true);
            if (do3DSRRF) queue.putReadBuffer(clBufferGz, true);

            bufferGx = clBufferGx.getBuffer();
            bufferGy = clBufferGy.getBuffer();
            if(do3DSRRF) bufferGz = clBufferGz.getBuffer();
        }
        else {
            imageWidth *= 2;
            imageHeight *= 2;
            queue.putReadBuffer(clBufferGxInt, true);
            queue.putReadBuffer(clBufferGyInt, true);
            if (do3DSRRF) queue.putReadBuffer(clBufferGzInt, true);

            bufferGx = clBufferGxInt.getBuffer();
            bufferGy = clBufferGyInt.getBuffer();
            if (do3DSRRF) bufferGz = clBufferGzInt.getBuffer();
        }

        ImageStack imsGradient = new ImageStack(imageWidth, imageHeight);
        if(do3DSRRF){
            // Load data
            for (int i = 0; i < nPlanesThreeD; i++) {

                float[] dataGx = new float[imageWidth * imageHeight];
                float[] dataGy = new float[imageWidth * imageHeight];
                float[] dataGz = new float[imageWidth * imageHeight];

                for (int n = 0; n < imageWidth * imageHeight; n++) {
                    dataGx[n] = bufferGx.get(n + i*imageWidth*imageHeight);
                    if (Float.isNaN(dataGx[n])) dataGx[n] = 0; // make sure we dont get any weirdness
                    dataGy[n] = bufferGy.get(n + i*imageWidth*imageHeight);
                    if (Float.isNaN(dataGy[n])) dataGy[n] = 0; // make sure we dont get any weirdness
                    dataGz[n] = bufferGz.get(n + i*imageWidth*imageHeight);
                    if (Float.isNaN(dataGz[n])) dataGz[n] = 0; // make sure we dont get any weirdness
                }

                imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGx));
                imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGy));
                imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGz));
            }

        }
        else { // while doing 2D
            // Load data
            float[] dataGx = new float[imageWidth * imageHeight];
            float[] dataGy = new float[imageWidth * imageHeight];

            for (int n = 0; n < imageWidth * imageHeight; n++) {
                dataGx[n] = bufferGx.get(n);
                if (Float.isNaN(dataGx[n])) dataGx[n] = 0; // make sure we dont get any weirdness
                dataGy[n] = bufferGy.get(n);
                if (Float.isNaN(dataGy[n])) dataGy[n] = 0; // make sure we dont get any weirdness
            }

            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGx));
            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGy));
        }

        return imsGradient;
    }


    // --- Release GPU context ---
    public void release() {
        //context.release();
        while (!context.isReleased()){
            IJ.log("-------------");
            IJ.log("Releasing context...");
            context.release();
        }
    }


    // --- Reset the SRRF frame counter ---
    public void resetFramePosition() {
        int id = prof.startTimer();
        kernelResetFramePosition.setArg(0, clBufferCurrentFrame); // make sure type is the same !!
        queue.put1DRangeKernel(kernelResetFramePosition, 0, 1, 0);
        prof.recordTime("Reset SRRF frame counter", prof.endTimer(id));
    }

    // --- Reshape the raw input data for 3D-SRRF ---
    public ImageStack reshapeImageStack3D(ImageStack ims) {

        int nPixels = ims.getWidth()*ims.getHeight()*ims.getSize();
        assert (widthS*heightS*nPlanesThreeD*ims.getSize() == nPixels);
        ImageStack reshapedIms = new ImageStack(widthS*heightS, nPlanesThreeD, ims.getSize());

        int z, xL, yL, xG, yG;
        for (int f = 0; f < ims.getSize(); f++) {

            FloatProcessor fp = ims.getProcessor(f+1).convertToFloatProcessor();
            float[] pixelsReshaped = new float[widthS*heightS*nPlanesThreeD];

            for (int i = 0; i < widthS*heightS*nPlanesThreeD; i++) {

                z = (int) chosenROIsLocationsThreeD[i/(widthS*heightS)]; // reorder the ROI to be in the right order
                yL = (i - z*widthS*heightS)/widthS; // local coordinates
                xL = i - z*widthS*heightS - yL*widthS;

                xG = xL + (z - 3*(z/3))*widthS; // global coordinate
                yG = yL + (z/3)*heightS;
                pixelsReshaped[i] = fp.getf(xG, yG)/ (float) intCoeffsThreeD[i/(widthS*heightS)]; // rescale the data according to the intensity factors

            }
            FloatProcessor fpOut = new FloatProcessor(widthS*heightS, nPlanesThreeD, pixelsReshaped);
            reshapedIms.setProcessor(fpOut, f+1);
        }

        // Checking the shape of the data if necessary
//        ImagePlus impReshaped = new ImagePlus("Reshaped data", reshapedIms);
//        impReshaped.show();

        return reshapedIms;
    }

}
