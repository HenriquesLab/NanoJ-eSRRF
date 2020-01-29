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

public class LiveSRRF_CL {

    // Basic formats
    private int width,
            height,
            widthM,
            heightM,
            nFrameForSRRF,
            blockLength,
            magnification,
            nFramesOnGPU,
            nFrameToLoad;

    private float vxy_offset;
    private int vxy_ArrayShift;

    public final int nReconstructions = 3; // Currently AVG, VAR (2nd order SOFI Tau=0) and 2nd order cumulants Tau=1
    public final String[] reconNames = new String[]{"AVG", "VAR", "TAC2"};

    public ImageStack imsSRRF;

    // Advanced formats
    private final NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programLiveSRRF;
    static private CLKernel kernelCalculateGradient,
            kernelIncrementFramePosition,
            kernelResetFramePosition,
            kernelCalculateSRRF,
            kernelCalculateVar;

    static private CLPlatform clPlatformMaxFlop;
    static private CLDevice clDeviceMaxFlop;
    static public CLDevice[] allCLdevices;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferDriftXY,
            clBufferPreviousFrame,
            clBufferOut;

    private CLBuffer<IntBuffer>
            clBufferCurrentFrame;


    // --- Constructor ---
    public LiveSRRF_CL() {
        // Nothing to see here. Keep calm and carry on. -----> -------->
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

    // ---------------------------------- Initialization method ----------------------------------
    public void initialise(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF, int blockLength, CLDevice chosenDevice, boolean intWeighting, String thisGradientChoice) {

        this.width = width;
        this.height = height;
        this.heightM = height * magnification;
        this.widthM = width * magnification;
        this.nFrameForSRRF = nFrameForSRRF;
        this.blockLength = blockLength;
        this.magnification = magnification;
        this.nFramesOnGPU = nFramesOnGPU;

        if (chosenDevice == null) {
//            IJ.log("Looking for the fastest device...");
            System.out.println("Using the fastest device...");
            context = CLContext.create(clPlatformMaxFlop);
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

        clBufferPx = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_ONLY);
        clBufferDriftXY = context.createFloatBuffer(2 * nFrameForSRRF, READ_ONLY);
        clBufferGx = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_WRITE);
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * widthM * heightM, READ_WRITE);
        clBufferPreviousFrame = context.createFloatBuffer(widthM * heightM, READ_WRITE);
        clBufferOut = context.createFloatBuffer((nReconstructions + 1) * widthM * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);
        // Current frame is a 2 element Int buffer:
        // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
        // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every tun of the method calculateSRRF (within the gradient calculation))

        switch (thisGradientChoice) {
            case "RobX":
                vxy_offset = 0.0f;
                vxy_ArrayShift = 0;
//            IJ.log("Using RobX");
                break;
            case "3pPlus":
            case "5pPlus":
                vxy_offset = 0.5f;
                vxy_ArrayShift = 0;
//            IJ.log("Using 3pPlus");
                break;
            case "3pX":
                vxy_offset = 0.5f;
                vxy_ArrayShift = 0;
//            IJ.log("Using 3pX");
                break;
        }

        // Create the program
        float sigma = fwhm / 2.354f;
        float analysisRadius = (float) ((int) (2 * sigma)) + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I

        String programString = getResourceAsString(LiveSRRF_CL.class, "liveSRRF_noInterpolation.cl");

        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
        programString = replaceFirst(programString, "$TWOSIGSQUARE$", "" + 2 * sigma*magnification * sigma*magnification);
        programString = replaceFirst(programString, "$TWOSIGpONE$", "" + (2 * sigma + 1)*magnification);

        programString = replaceFirst(programString, "$RADIUS$", "" + magnification*analysisRadius);
        programString = replaceFirst(programString, "$WIDTH$", "" + widthM);
        programString = replaceFirst(programString, "$HEIGHT$", "" + heightM);
        programString = replaceFirst(programString, "$WH$", "" + (widthM * heightM));
        programString = replaceFirst(programString, "$VXY_OFFSET$", "" + vxy_offset);
        programString = replaceFirst(programString, "$VXY_ARRAYSHIFT$", "" + vxy_ArrayShift);
        programString = replaceFirst(programString, "$NFRAMEFORSRRF$", "" + nFrameForSRRF);

        if (intWeighting) programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 1);
        else programString = replaceFirst(programString, "$INTWEIGHTING$", "" + 0);

        programLiveSRRF = context.createProgram(programString).build();
//        IJ.log("Program executable? "+programLiveSRRF.isExecutable());
//        IJ.log("------------------------------------");
//        IJ.log(programLiveSRRF.getBuildLog());
//        IJ.log("------------------------------------");

        switch (thisGradientChoice) {
            case "RobX":
                kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradientRobX");
                break;
            case "3pPlus":
                kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient3pPlus");
                break;
            case "3pX":
                kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient3pX");
                break;
            case "5pPlus":
                kernelCalculateGradient = programLiveSRRF.createCLKernel(("calculateGradient5pPlus"));
                break;
        }

        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateRadialGradientConvergence");
        kernelIncrementFramePosition = programLiveSRRF.createCLKernel("kernelIncrementFramePosition");
        kernelResetFramePosition = programLiveSRRF.createCLKernel("kernelResetFramePosition");
        kernelCalculateVar = programLiveSRRF.createCLKernel("kernelCalculateVar");

        int argn;
        argn = 0;
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferPreviousFrame);
        kernelCalculateSRRF.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferDriftXY); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelIncrementFramePosition.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelCalculateVar.setArg(argn++, clBufferOut);

        queue = chosenDevice.createCommandQueue();

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferDriftXY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferPreviousFrame.getCLSize() +
                        clBufferOut.getCLSize() +
                        clBufferCurrentFrame.getCLSize()
        )
                / 1000000d + "MB"); // conversion in base 10
    }


    // --- Load Drift array on GPU ---
    public void loadDriftXYGPUbuffer(float[][] driftXY) {

//        System.out.println("loading drift buffers");
        float[] driftXYarray = new float[2 * nFrameForSRRF];
        for (int i = 0; i < nFrameForSRRF; i++) {
            driftXYarray[i] = magnification*driftXY[i][0]; // loaded as drift in pixels in the magnified pixel space
            driftXYarray[i + nFrameForSRRF] = magnification*driftXY[i][1];
        }

        int id = prof.startTimer();
        fillBuffer(clBufferDriftXY, driftXYarray);
        queue.putWriteBuffer(clBufferDriftXY, false);
        prof.recordTime("Uploading drift array to GPU", prof.endTimer(id));
    }


    // --- Calculate SRRF images ---
    public void prepareDataSRRF(ImageStack imsRawData) {

        assert (imsRawData.getWidth() == width && imsRawData.getHeight() == height);
        nFrameToLoad = imsRawData.getSize();

        int id;
        if (magnification > 1) {
            id = prof.startTimer();
            ImageStack imsRawDataInt = fhtSpaceInterpolation(imsRawData, magnification, true); // enables mirror padding by default
            prof.recordTime("Interpolating raw data", prof.endTimer(id));

//        IJ.log("Uploading raw data to GPU...");
            id = prof.startTimer();
            fillBuffer(clBufferPx, imsRawDataInt);
            queue.putWriteBuffer(clBufferPx, false);
            prof.recordTime("Uploading data to GPU", prof.endTimer(id));
        } else {

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
    }

    // --- Calculate SRRF images ---
    public synchronized boolean calculateSRRF() { // returns boolean describing whether it was cancelled by user or not

        int workSize;
        int nBlocks = widthM * heightM / blockLength + ((widthM * heightM % blockLength == 0) ? 0 : 1);
        queue.finish(); // Make sure everything is done

        int id;
        for (int f = 0; f < nFrameToLoad; f++) {
//            System.out.println("Calculate SRRF frame #"+f);

            for (int nB = 0; nB < nBlocks; nB++) {
                workSize = min(blockLength, widthM * heightM - nB * blockLength);
                id = prof.startTimer();
                queue.put1DRangeKernel(kernelCalculateSRRF, nB * blockLength, workSize, 0); // TODO: compare 1D kernel vs ND kernel
                prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return true;
                }
            }

            // This kernel needs to be done outside of the previous kernel because of concomitent execution (you never know when each pixel is executed)
            id = prof.startTimer();
            queue.finish(); // Make sure everything is done
            queue.put1DRangeKernel(kernelIncrementFramePosition, 0, 2, 0); // this internally increment the frame position f
            prof.recordTime("Increment frame count", prof.endTimer(id));
        }

        return false;

    }

    // --- Read the output buffer ---
    public void readSRRFbuffer() {

//        System.out.println("Reading SRRF buffer");
        // Calculate the VARIANCE on the OutputArray on the GPU
        int id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put1DRangeKernel(kernelCalculateVar, 0, heightM*widthM,0);
        prof.recordTime("Calculate VAR image", prof.endTimer(id));

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferOut, true);
        FloatBuffer bufferSRRF = clBufferOut.getBuffer();

        imsSRRF = new ImageStack(widthM, heightM);

        // Load data back to CPU memory
        for (int i = 0; i < nReconstructions+1; i++) {
            float[] dataSRRF = new float[widthM * heightM];
            for (int n = 0; n < widthM * heightM; n++) {
                dataSRRF[n] = bufferSRRF.get(n + i*widthM * heightM);
                if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we don't get any weirdness
            }
            imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));
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

//        System.out.println("Reading Gx/Gy buffers");
        queue.putReadBuffer(clBufferGx, true);
        queue.putReadBuffer(clBufferGy, true);
//        System.out.println("Gx/Gy buffers were read successfully");

        bufferGx = clBufferGx.getBuffer();
        bufferGy = clBufferGy.getBuffer();

        ImageStack imsGradient = new ImageStack(imageWidth, imageHeight);
        // Load data
        for (int i = 0; i < nFramesOnGPU; i++) {

            float[] dataGx = new float[imageWidth * imageHeight];
            float[] dataGy = new float[imageWidth * imageHeight];

            for (int n = 0; n < imageWidth * imageHeight; n++) {
                dataGx[n] = bufferGx.get(n + i*imageWidth*imageHeight);
                if (Float.isNaN(dataGx[n])) dataGx[n] = 0; // make sure we dont get any weirdness
                dataGy[n] = bufferGy.get(n + i*imageWidth*imageHeight);
                if (Float.isNaN(dataGy[n])) dataGy[n] = 0; // make sure we dont get any weirdness
            }

            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGx));
            imsGradient.addSlice(new FloatProcessor(imageWidth, imageHeight, dataGy));
        }

        return imsGradient;
    }


    // --- Release GPU context ---
    public void release() {

        if (context != null) {
            while (!context.isReleased()) {
                IJ.log("-------------");
                IJ.log("Releasing context...");
                context.release();
            }
        }
    }

}
