package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.getResourceAsString;
import static nanoj.core2.NanoJCL.replaceFirst;
import static nanoj.liveSRRF.FourierInterpolation.fhtSpaceInterpolation;

public class liveSRRF_CL {

    // Basic formats
    private int width,
            height,
            widthM,
            heightM,
            nFrameForSRRF,
            blockLength,
            magnification,
            nFramesOnGPU;

    private final int GxGyMagnification = 1;
    //    private final float vxy_offset = 0.5f; // These are for the 2pI method
//    private final int vxy_ArrayShift = 1;
//    private final float vxy_offset = 0.0f; //These are for the RobX method
//    private final int vxy_ArrayShift = 0;
    private float vxy_offset; //These are for the RobX method
    private int vxy_ArrayShift;

    private final int nReconstructions = 2; // Currently only STD and AVG

    private boolean doMPmapCorrection;

    public ImageStack imsSRRF;

    private String thisGradientChoice;

    // Advanced formats
    private final NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programLiveSRRF;
    static private CLKernel kernelCalculateGradient,
    //                            kernelInterpolateGradient,
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
            clBufferGx, clBufferGy,
    //            clBufferGxInt, clBufferGyInt,
    clBufferShiftXY,
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
    public void initialise(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF, int blockLength, CLDevice chosenDevice, boolean intWeighting, boolean doMPmapCorrection, String thisGradientChoice) {

        this.width = width;
        this.height = height;
        this.heightM = height * magnification;
        this.widthM = width * magnification;
        this.nFrameForSRRF = nFrameForSRRF;
        this.blockLength = blockLength;
        this.magnification = magnification;
        this.doMPmapCorrection = doMPmapCorrection;
        this.thisGradientChoice = thisGradientChoice;
        this.nFramesOnGPU = nFramesOnGPU;

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


        System.out.println("using " + chosenDevice);
        //IJ.log("Using " + chosenDevice.getName());

//        clBufferPx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_ONLY);
        clBufferPx = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_ONLY);
        clBufferShiftXY = context.createFloatBuffer(2 * nFrameForSRRF, READ_ONLY);
//        clBufferGx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gx
//        clBufferGy = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gy

        clBufferGx = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_WRITE); // single frame Gx
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_WRITE); // single frame Gy


//        clBufferGxInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gx
//        clBufferGyInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gy
        clBufferOut = context.createFloatBuffer((nReconstructions + 1) * widthM * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);
        clBufferMPmap = context.createFloatBuffer(2 * magnification * magnification, READ_WRITE);

        // Current frame is a 2 element Int buffer:
        // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
        // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every tun of the method calculateSRRF (within the gradient calculation))

        if (thisGradientChoice.equals("RobX")) {
            vxy_offset = 0.0f;
            vxy_ArrayShift = 0;
            IJ.log("Using RobX");
        }
        else if (thisGradientChoice.equals("3pPlus")) {
            vxy_offset = 0.5f;
            vxy_ArrayShift = 0;
            IJ.log("Using 3pPlus");
        }
        else if (thisGradientChoice.equals("3pX")) {
            vxy_offset = 0.5f;
            vxy_ArrayShift = 0;
            IJ.log("Using 3pX");
        }

        // Create the program
        float sigma = fwhm / 2.354f;
        float radius = ((float) ((int) (GxGyMagnification * 2 * sigma))) / GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I
//        float radius = 2*((float) ((int) (GxGyMagnification * 2 * sigma))) / GxGyMagnification + 1;    // TODO: for testing


        String programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");
//        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
//        programString = replaceFirst(programString, "$FWHM$", "" + fwhm);
//        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
//        programString = replaceFirst(programString, "$GXGYMAGNIFICATION$", "" + GxGyMagnification);

//        programString = replaceFirst(programString, "$SIGMA$", "" + sigma);
//        programString = replaceFirst(programString, "$RADIUS$", "" + radius);
//        programString = replaceFirst(programString, "$WIDTH$", "" + width);
//        programString = replaceFirst(programString, "$HEIGHT$", "" + height);
//        programString = replaceFirst(programString, "$WH$", "" + (width * height));
//        programString = replaceFirst(programString, "$WM$", "" + (width * magnification));
//        programString = replaceFirst(programString, "$HM$", "" + (height * magnification));
//        programString = replaceFirst(programString, "$WHM$", "" + (width * height * magnification * magnification));
//        programString = replaceFirst(programString, "$WINT$", "" + (GxGyMagnification * width));
//        programString = replaceFirst(programString, "$HINT$", "" + (GxGyMagnification * height));
//        programString = replaceFirst(programString, "$VXY_OFFSET$", "" + vxy_offset);
//        programString = replaceFirst(programString, "$VXY_ARRAYSHIFT$", "" + vxy_ArrayShift);
//        programString = replaceFirst(programString, "$NFRAMEFORSRRF$", "" + nFrameForSRRF);

        // This lot of definition works for no interpolation on GPU only (GxGyMagnification = 1 and width/height widthM/heightM are the same from the GPU's point of view)
        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + 1);
        programString = replaceFirst(programString, "$FWHM$", "" + magnification*fwhm);
        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
        programString = replaceFirst(programString, "$GXGYMAGNIFICATION$", "" + 1);

        programString = replaceFirst(programString, "$SIGMA$", "" + magnification*sigma);
        programString = replaceFirst(programString, "$RADIUS$", "" + magnification*radius);
        programString = replaceFirst(programString, "$WIDTH$", "" + widthM);
        programString = replaceFirst(programString, "$HEIGHT$", "" + heightM);
        programString = replaceFirst(programString, "$WH$", "" + (widthM * heightM));
        programString = replaceFirst(programString, "$WM$", "" + widthM);
        programString = replaceFirst(programString, "$HM$", "" + heightM);
        programString = replaceFirst(programString, "$WHM$", "" + (widthM * heightM));
        programString = replaceFirst(programString, "$WINT$", "" + widthM);
        programString = replaceFirst(programString, "$HINT$", "" + heightM);
        programString = replaceFirst(programString, "$VXY_OFFSET$", "" + vxy_offset);
        programString = replaceFirst(programString, "$VXY_ARRAYSHIFT$", "" + vxy_ArrayShift);
        programString = replaceFirst(programString, "$NFRAMEFORSRRF$", "" + nFrameForSRRF);

        if (intWeighting) programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 1);
        else programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 0);

        programLiveSRRF = context.createProgram(programString).build();

//        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");

        if (thisGradientChoice.equals("RobX")) kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradientRobX");
        else if (thisGradientChoice.equals("3pPlus")) kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient3pPlus");
        else if (thisGradientChoice.equals("3pX")) kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient3pX");

//        kernelInterpolateGradient = programLiveSRRF.createCLKernel("calculateGradientInterpolation");
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
        kernelCalculateGradient.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

//        argn = 0;
//        kernelInterpolateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
//        kernelInterpolateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
//        kernelInterpolateGradient.setArg(argn++, clBufferGxInt); // make sure type is the same !!
//        kernelInterpolateGradient.setArg(argn++, clBufferGyInt); // make sure type is the same !!

        argn = 0;
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGy); // make sure type is the same !!
//        kernelCalculateSRRF.setArg(argn++, clBufferGxInt); // make sure type is the same !!
//        kernelCalculateSRRF.setArg(argn++, clBufferGyInt); // make sure type is the same !!
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
//                        clBufferGxInt.getCLSize() +
//                        clBufferGyInt.getCLSize() +
                        clBufferOut.getCLSize() +
                        clBufferCurrentFrame.getCLSize() +
                        clBufferMPmap.getCLSize()
        )
                / 1000000d + "MB");

    }


    // --- Load Shift array on GPU ---
    public void loadShiftXYGPUbuffer(float[] shiftX, float[] shiftY) {

//        System.out.println("loading shift buffers");

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

//        System.out.println("Calculate SRRF entry");

        assert (imsRawData.getWidth() == width && imsRawData.getHeight() == height);
        int nFrameToLoad = imsRawData.getSize();

        int id;
        if (magnification > 1) {
            id = prof.startTimer();
            ImageStack imsRawDataInt = fhtSpaceInterpolation(imsRawData, magnification, true); // enables mirror padding by default
            prof.recordTime("Interpolating raw data", prof.endTimer(id));

//        IJ.log("Uploading raw data to GPU...");
            id = prof.startTimer();
//        fillBuffer(clBufferPx, imsRawData);
            fillBuffer(clBufferPx, imsRawDataInt);
            queue.putWriteBuffer(clBufferPx, false);
            prof.recordTime("Uploading data to GPU", prof.endTimer(id));
        }
        else {

//        IJ.log("Uploading raw data to GPU...");
            id = prof.startTimer();
            fillBuffer(clBufferPx, imsRawData);
            queue.putWriteBuffer(clBufferPx, false);
            prof.recordTime("Uploading data to GPU", prof.endTimer(id));
        }

        // Make kernelCalculateGradient assignment (this kernel also resets the local GPU load frame counter)
//        IJ.log("Calculating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put3DRangeKernel(kernelCalculateGradient, 0, 0, 0, widthM, heightM, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

////        IJ.log("Interpolating gradient...");
//        id = prof.startTimer();
//        queue.finish(); // Make sure everything is done
//        queue.put3DRangeKernel(kernelInterpolateGradient, 0, 0, 0, GxGyMagnification * width, GxGyMagnification * height, nFrameToLoad, 0, 0, 0);
//        prof.recordTime("kernelInterpolateGradient", prof.endTimer(id));

        // Make kernelCalculateSRRF assignment
//        IJ.log("Calculating SRRF...");

        int workSize;
        int nBlocks = widthM * heightM / blockLength + ((widthM * heightM % blockLength == 0) ? 0 : 1);

        for (int f = 0; f < nFrameToLoad; f++) {
//            System.out.println("Calculate SRRF frame #"+f);

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

//        System.out.println("Reading SRRF buffer");

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

    public void checkProfiler(){
        IJ.log("-------------------------------------");
        IJ.log(prof.report());
    }

    // --- Read the gradient buffers ---
    public ImageStack readGradientBuffers() {

        queue.finish(); // Make sure everything is done

        int imageWidth = widthM;
        int imageHeight = heightM;

        FloatBuffer bufferGx;
        FloatBuffer bufferGy;
//        FloatBuffer bufferGz = null;

//        if (!interpolated){
        System.out.println("Reading Gx/Gy buffers");
        queue.putReadBuffer(clBufferGx, true);
        queue.putReadBuffer(clBufferGy, true);
        System.out.println("Gx/Gy buffers were read successfully");


        bufferGx = clBufferGx.getBuffer();
        bufferGy = clBufferGy.getBuffer();
//        }
//        else {
//            imageWidth *= gradientMag;
//            imageHeight *= gradientMag;
//
//            queue.putReadBuffer(clBufferGxInt, true);
//            queue.putReadBuffer(clBufferGyInt, true);
//            if (do3DSRRF) queue.putReadBuffer(clBufferGzInt, true);
//
//            bufferGx = clBufferGxInt.getBuffer();
//            bufferGy = clBufferGyInt.getBuffer();
//            if (do3DSRRF) bufferGz = clBufferGzInt.getBuffer();
//        }

        ImageStack imsGradient = new ImageStack(imageWidth, imageHeight);
        // Load data // TODO: this currently only outputs the first time frame but all the planes
        for (int i = 0; i < nFramesOnGPU; i++) {

            float[] dataGx = new float[imageWidth * imageHeight];
            float[] dataGy = new float[imageWidth * imageHeight];

            for (int n = 0; n < imageWidth * imageHeight; n++) {
                dataGx[n] = bufferGx.get(n + i*imageWidth*imageHeight);
                if (Float.isNaN(dataGx[n])) dataGx[n] = 0; // make sure we dont get any weirdness
                dataGy[n] = bufferGy.get(n + i*imageWidth*imageHeight);
                if (Float.isNaN(dataGy[n])) dataGy[n] = 0; // make sure we dont get any weirdness
//                if (do3DSRRF) {
//                    dataGz[n] = bufferGz.get(n + i * imageWidth * imageHeight);
//                    if (Float.isNaN(dataGz[n])) dataGz[n] = 0; // make sure we dont get any weirdness
//                }
            }

            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGx));
            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGy));
//            if (do3DSRRF) imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGz));
        }

        return imsGradient;
    }

}
