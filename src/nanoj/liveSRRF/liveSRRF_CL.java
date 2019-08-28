package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import nanoj.core.java.array.ArrayMath;
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
import static nanoj.liveSRRF.LinearRegressions.linearRegressionLeastSquareNoOffset;

public class liveSRRF_CL {

    // Basic formats
    private int singleFrameSize,
            widthM,
            heightM,
            nFrameForSRRF,
            blockLength,
            magnification,
            nSplits = 3, // TODO: currently hard-coded
            nPlanes;

    public int widthS, heightS;

    private double[] intCoeffs3D,
            chosenROIsLocations3D;

    public final int gradientMag = 2;
    public final int nReconstructions = 2; // Currently only STD and AVG

    private final float vxy_offset = 0.5f;
    private final int vxy_ArrayShift = 1;

    private float dimensionRatioXYvsZ;

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
            clBufferDriftXY,
            clBufferShiftXY3D, // TODO: add the shift for 3D-SRRF!
            clBufferOut,
            clBufferMPmap;

    private CLBuffer<IntBuffer>
            clBufferCurrentFrame;

    // ---------------------------------- Constructor ----------------------------------
    public liveSRRF_CL() {
        // Nothing to see here. Keep calm and carry on.
    }

    // ---------------------------------- Check devices ----------------------------------
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


    // ---------------------------------- Initialization method ----------------------------------
    public void initialise(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF, int blockLength, CLDevice chosenDevice, boolean intWeighting, boolean doMPmapCorrection, String calibTablePath3DSRRF, float axialOffset, float pixelSize) {

//        this.width = width;
//        this.height = height;
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
        else {
            context = CLContext.create(chosenDevice.getPlatform());
            IJ.log("Working on platform: "+chosenDevice.getPlatform().getName());
            CLDevice[] allCLdevicesOnThisPlatform = context.getDevices();
            int i = 0;
            while (!allCLdevicesOnThisPlatform[i].getName().equals(chosenDevice.getName())) {
                i++;
            }
            chosenDevice = allCLdevicesOnThisPlatform[i];
        }

        // 3D-SRRF initialisation
        do3DSRRF = (calibTablePath3DSRRF != null);
        if (do3DSRRF) IJ.log("3D-SRRRRRRRRRF!!!");
//        System.out.print("voila!!");

        double[] shiftX3D = null;
        double[] shiftY3D = null;
        double[] theta3D = null; // TODO: add the theta correction to the recon

        if (do3DSRRF) {
            double[] axialPositions3D = null;
            double[] nominalPositions3D = null;

            Map<String, double[]> calibTable;
            try {
                calibTable = new LoadNanoJTable(calibTablePath3DSRRF).getData();
                shiftX3D = calibTable.get("X-shift (pixels)");
                shiftY3D = calibTable.get("Y-shift (pixels)");
                theta3D = calibTable.get("Theta (degrees)");
                axialPositions3D = calibTable.get("Axial positions");
                nominalPositions3D = calibTable.get("Nominal positions");
                chosenROIsLocations3D = calibTable.get("ROI #");
                intCoeffs3D = calibTable.get("Intensity scaling");
                ResultsTable rt = dataMapToResultsTable(calibTable);
                rt.show("Calibration-Table");
            } catch (IOException e) {
                IJ.log("Catching exception...");
                e.printStackTrace();
            }
            nPlanes = shiftX3D.length;
            widthS = width / nSplits;
            heightS = height / nSplits;

            if (ArrayMath.sum(axialPositions3D) != 0)
                axialOffset = (float) linearRegressionLeastSquareNoOffset(nominalPositions3D, axialPositions3D);
            dimensionRatioXYvsZ = axialOffset / (1000 * pixelSize); // TODO: this assumes that pixelSize is in um!
            IJ.log("3D axial offset: " + axialOffset + " nm");
            IJ.log("Dimension ratio (Z/XY) " + dimensionRatioXYvsZ);

        } else { // if doing 2D
            nPlanes = 1;
            widthS = width;
            heightS = height;
        }

        singleFrameSize = nPlanes * widthS * heightS; // nPlanes = 1 in 2D case
        this.widthM = widthS * magnification;
        this.heightM = heightS * magnification;
        IJ.log("WidthM/HeightM: "+widthM+"/"+heightM);

        System.out.println("using " + chosenDevice);
        //IJ.log("Using " + chosenDevice.getName());

        // initialise buffers
        clBufferPx = context.createFloatBuffer(nFramesOnGPU * singleFrameSize, READ_ONLY);
        clBufferDriftXY = context.createFloatBuffer(2 * nFrameForSRRF, READ_ONLY);
        if (do3DSRRF) clBufferShiftXY3D = context.createFloatBuffer(2 * nPlanes, READ_ONLY);
        clBufferGx = context.createFloatBuffer(nFramesOnGPU * singleFrameSize, READ_WRITE); // all frames Gx
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * singleFrameSize, READ_WRITE); // all frames Gy
        clBufferGxInt = context.createFloatBuffer(nFramesOnGPU * gradientMag * gradientMag * singleFrameSize, READ_WRITE); // single frame Gx
        clBufferGyInt = context.createFloatBuffer(nFramesOnGPU * gradientMag * gradientMag * singleFrameSize, READ_WRITE); // single frame Gy

        // initialise buffers for 3D-SRRF
        if (do3DSRRF) {
            clBufferShiftXY3D = context.createFloatBuffer(2 * nPlanes, READ_ONLY);
            // TODO: add some buffers for the rotation
            clBufferGz = context.createFloatBuffer(nFramesOnGPU * singleFrameSize, READ_WRITE);
            clBufferGzInt = context.createFloatBuffer(nFramesOnGPU * gradientMag * gradientMag * singleFrameSize, READ_WRITE);
            // TODO: this assumes only interpolation in XY, not in Z yet
        }

        if (do3DSRRF)
            clBufferOut = context.createFloatBuffer((nReconstructions + 1) * singleFrameSize * magnification * magnification * magnification, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        else
            clBufferOut = context.createFloatBuffer((nReconstructions + 1) * singleFrameSize * magnification * magnification, WRITE_ONLY); // single frame cumulative AVG projection of RGC

        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);
        // Current frame is a 2 element Int buffer:
        // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
        // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every tun of the method calculateSRRF (within the gradient calculation))
        clBufferMPmap = context.createFloatBuffer(2 * magnification * magnification, READ_WRITE);

        // Create the program
        float sigma = fwhm / 2.354f;
        float radius = ((float) ((int) (gradientMag * 2 * sigma))) / gradientMag + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I

        String programString;
        if (do3DSRRF) programString = getResourceAsString(liveSRRF_CL.class, "live3DSRRF.cl");
        else programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");

        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
        programString = replaceFirst(programString, "$FWHM$", "" + fwhm);
        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
        programString = replaceFirst(programString, "$GXGYMAGNIFICATION$", "" + gradientMag);

        programString = replaceFirst(programString, "$SIGMA$", "" + sigma);
        programString = replaceFirst(programString, "$SIGMA22$", "" + 2 * sigma * sigma);
        programString = replaceFirst(programString, "$RADIUS$", "" + radius);
        programString = replaceFirst(programString, "$WIDTH$", "" + widthS);
        programString = replaceFirst(programString, "$HEIGHT$", "" + heightS);
        programString = replaceFirst(programString, "$WH$", "" + (widthS * heightS));
        programString = replaceFirst(programString, "$WM$", "" + (widthS * magnification));
        programString = replaceFirst(programString, "$HM$", "" + (heightS * magnification));
        programString = replaceFirst(programString, "$WHM$", "" + (widthS * heightS * magnification * magnification));
        programString = replaceFirst(programString, "$WINT$", "" + (gradientMag * widthS));
        programString = replaceFirst(programString, "$HINT$", "" + (gradientMag * heightS));
        programString = replaceFirst(programString, "$WHINT$", "" + (gradientMag * widthS * gradientMag * heightS));
        programString = replaceFirst(programString, "$VXY_OFFSET$", "" + vxy_offset);
        programString = replaceFirst(programString, "$VXY_ARRAYSHIFT$", "" + vxy_ArrayShift);
        programString = replaceFirst(programString, "$NFRAMEFORSRRF$", "" + nFrameForSRRF);

        if (do3DSRRF) {
            programString = replaceFirst(programString, "$NPLANES$", "" + nPlanes);
            programString = replaceFirst(programString, "$DIMRATIO$", "" + dimensionRatioXYvsZ);
            programString = replaceFirst(programString, "$WHD$", "" + widthS * heightS * nPlanes);
            programString = replaceFirst(programString, "$WHDM$", "" + widthS * heightS * nPlanes * magnification * magnification * magnification);
        }

        if (intWeighting) programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 1);
        else programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 0);

        programLiveSRRF = context.createProgram(programString).build();
        IJ.log("Program executable? "+programLiveSRRF.isExecutable());
        IJ.log(programLiveSRRF.getBuildLog());

        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");
        kernelInterpolateGradient = programLiveSRRF.createCLKernel("calculateGradientInterpolation");
        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateRadialGradientConvergence");
        kernelIncrementFramePosition = programLiveSRRF.createCLKernel("kernelIncrementFramePosition");
        kernelResetFramePosition = programLiveSRRF.createCLKernel("kernelResetFramePosition");
        kernelCalculateMPmap = programLiveSRRF.createCLKernel("kernelCalculateMPmap");
        kernelCalculateStd = programLiveSRRF.createCLKernel("kernelCalculateStd");
        kernelCorrectMPmap = programLiveSRRF.createCLKernel("kernelCorrectMPmap");
        IJ.log("Program executable? "+programLiveSRRF.isExecutable());

        int argn;
        argn = 0;
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        if (do3DSRRF) kernelCalculateGradient.setArg(argn++, clBufferGz); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelInterpolateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        if (do3DSRRF) kernelInterpolateGradient.setArg(argn++, clBufferGz); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        if (do3DSRRF) kernelInterpolateGradient.setArg(argn++, clBufferGzInt); // make sure type is the same !!

        argn = 0;
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        if (do3DSRRF) kernelCalculateSRRF.setArg(argn++, clBufferGzInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferDriftXY); // make sure type is the same !!
        if (do3DSRRF) kernelCalculateSRRF.setArg(argn++, clBufferShiftXY3D); // make sure type is the same !!
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
                        clBufferDriftXY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferGxInt.getCLSize() +
                        clBufferGyInt.getCLSize() +
                        clBufferOut.getCLSize() +
                        clBufferCurrentFrame.getCLSize() +
                        clBufferMPmap.getCLSize() // TODO: add the stuff from 3D
        )
                / 1000000d + "MB");


        if (do3DSRRF) {
            // Load the buffer for the XY shift for 3D
            float[] shiftXY3D = new float[2 * nPlanes];
            for (int i = 0; i < 2 * nPlanes; i++) {
                if (i < nPlanes) shiftXY3D[i] = (float) shiftX3D[i];
                if (i >= nPlanes) shiftXY3D[i] = (float) shiftY3D[i-nPlanes];
            }

            int id = prof.startTimer();
            fillBuffer(clBufferShiftXY3D, shiftXY3D);
            queue.putWriteBuffer(clBufferShiftXY3D, false);
            prof.recordTime("Uploading ShiftXY3D array to GPU", prof.endTimer(id));
        }

        IJ.log("Program executable? "+programLiveSRRF.isExecutable());

    }


    // --- Load Drift array on GPU ---
    public void loadDriftXYGPUbuffer(float[] driftX, float[] driftY) {

        System.out.println("TAG1");
        float[] driftXY = new float[2 * nFrameForSRRF]; // TODO: change to do manually
        System.arraycopy(driftX, 0, driftXY, 0, nFrameForSRRF);
        System.arraycopy(driftY, 0, driftXY, nFrameForSRRF, nFrameForSRRF);

        int id = prof.startTimer();
        fillBuffer(clBufferDriftXY, driftXY);
        queue.putWriteBuffer(clBufferDriftXY, false);
        prof.recordTime("Uploading drift array to GPU", prof.endTimer(id));
    }


    // --- Calculate SRRF images ---
    public synchronized boolean calculateSRRF(ImageStack imsRawData) { // returns boolean describing whether it was cancelled by user or not

        //assert (imsRawData.getWidth() == width && imsRawData.getHeight() == height); // TODO: put it back?
        int nFrameToLoad = imsRawData.getSize();

        //        IJ.log("Uploading raw data to GPU...");
        int id = prof.startTimer();
        if(do3DSRRF){
            ImageStack imsData3D = reshapeImageStack3D(imsRawData);
            fillBuffer(clBufferPx, imsData3D);
        }
        else fillBuffer(clBufferPx, imsRawData);
        System.out.println("TAG2");

        queue.putWriteBuffer(clBufferPx, false);
        prof.recordTime("Uploading data to GPU", prof.endTimer(id));

        // Make kernelCalculateGradient assignment (this kernel also resets the local GPU load frame counter)
//        IJ.log("Calculating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        if(do3DSRRF) queue.put1DRangeKernel(kernelCalculateGradient, 0,  singleFrameSize * nFrameToLoad, 0);
        else         queue.put3DRangeKernel(kernelCalculateGradient, 0, 0, 0, widthS, heightS, nFrameToLoad, 0, 0, 0);
        System.out.println("TAG3");

        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

//        IJ.log("Interpolating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        if(do3DSRRF) queue.put1DRangeKernel(kernelInterpolateGradient, 0,  gradientMag * gradientMag * singleFrameSize * nFrameToLoad, 0);
        else         queue.put3DRangeKernel(kernelInterpolateGradient, 0, 0, 0, gradientMag * widthS, gradientMag * heightS, nFrameToLoad, 0, 0, 0);

        prof.recordTime("kernelInterpolateGradient", prof.endTimer(id));
        System.out.println("TAG4");

        // Make kernelCalculateSRRF assignment, in blocks as defined by user (blocksize from GUI)
//        IJ.log("Calculating SRRF...");

        int workSize;

        // Calculate the number of blocks necessary to compute the whole thing
        int nBlocks;
        if (do3DSRRF) nBlocks = nPlanes * magnification * widthM * heightM / blockLength + ((nPlanes * magnification * widthM * heightM % blockLength == 0) ? 0 : 1);
        else nBlocks = widthM * heightM / blockLength + ((widthM * heightM % blockLength == 0) ? 0 : 1);
        System.out.println("TAG5");

        for (int f = 0; f < nFrameToLoad; f++) {
            for (int nB = 0; nB < nBlocks; nB++) {
                // Calculate the worksize, always the blockLength until there's fewer pixels left than blocksize
                // TODO: optimise to get an exact integer of total number of pixels and minimize the number of calls?
                if (do3DSRRF) workSize = min(blockLength, nPlanes * magnification * widthM * heightM - nB * blockLength);
                else workSize = min(blockLength, widthM * heightM - nB * blockLength);

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
            queue.put1DRangeKernel(kernelIncrementFramePosition, 0, 2, 0); // this internally increment the frame position f
            prof.recordTime("Increment frame count", prof.endTimer(id));
        }

        System.out.println("TAG6");

        return false;

    }

    // --- Read the output buffer ---
    public void readSRRFbuffer() {

        // Calculate the STD on the OutputArray on the GPU
        int id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put1DRangeKernel(kernelCalculateStd, 0, singleFrameSize*magnification*magnification,0);
        prof.recordTime("Calculate STD image", prof.endTimer(id));

        if (doMPmapCorrection) { // TODO: is this going to be a problem for 3D? probably yes
            id = prof.startTimer();
            queue.put1DRangeKernel(kernelCalculateMPmap, 0, 2 * magnification * magnification, 0);
            prof.recordTime("Calculate MP map", prof.endTimer(id));

            id = prof.startTimer();
            queue.put1DRangeKernel(kernelCorrectMPmap, 0, 2 * singleFrameSize*magnification*magnification, 0);
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

        int imageWidth = widthS;
        int imageHeight = heightS;

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
            for (int i = 0; i < nPlanes; i++) {

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
        assert (widthS*heightS* nPlanes *ims.getSize() == nPixels);
        ImageStack reshapedIms = new ImageStack(widthS*heightS, nPlanes, ims.getSize());

        int z, xL, yL, xG, yG, offset;
        for (int f = 0; f < ims.getSize(); f++) {

            FloatProcessor fp = ims.getProcessor(f+1).convertToFloatProcessor();
            float[] pixelsReshaped = new float[widthS*heightS* nPlanes];

            for (int i = 0; i < widthS*heightS* nPlanes; i++) {

                z = i/(widthS*heightS); // reorder the ROI to be in the right order
//                IJ.log("Reshape z: "+z);
                yL = (i - z*widthS*heightS)/widthS; // local coordinates
                xL = i - z*widthS*heightS - yL*widthS;
                xG = xL + (z - 3*(z/3))*widthS; // global coordinate
                yG = yL + (z/3)*heightS;

                offset = (int) chosenROIsLocations3D[z]*widthS*heightS + yL*widthS + xL;
//                IJ.log("Reshape offset: "+offset);
                pixelsReshaped[offset] = fp.getf(xG, yG) / (float) intCoeffs3D[z]; // rescale the data according to the intensity factors
            }
            FloatProcessor fpOut = new FloatProcessor(widthS*heightS, nPlanes, pixelsReshaped);
            reshapedIms.setProcessor(fpOut, f+1);
        }

        // Checking the shape of the data if necessary
//        ImagePlus impReshaped = new ImagePlus("Reshaped data", reshapedIms);
//        impReshaped.show();

        return reshapedIms;
    }

}
