package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;
import org.python.modules.math;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.min;
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
            magnification;

    private final int GxGyMagnification = 2;
    private final float vxy_offset = 0.5f;
    private final int vxy_ArrayShift = 1;

    private final int nReconstructions = 2; // Currently only STD and AVG

    private final boolean DEBUG = false;

    public ImageStack imsSRRF;

    // Advanced formats
    private NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programLiveSRRF;
    static private CLKernel kernelCalculateGradient,
                            kernelInterpolateGradient,
                            kernelIncrementFramePosition,
                            kernelResetFramePosition,
                            kernelCalculateSRRF,
                            kernelCalculateMPmap,
                            kernelCalculateStd;

    static private CLPlatform clPlatformMaxFlop;
    static private CLDevice clDeviceMaxFlop;
    static public CLDevice[] allCLdevices;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferGxInt, clBufferGyInt,
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
    public void initialise(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF, int blockLength, CLDevice chosenDevice, boolean intWeighting) {

        this.width = width;
        this.height = height;
        this.heightM = height * magnification;
        this.widthM = width * magnification;
        this.nFrameForSRRF = nFrameForSRRF;
        this.blockLength = blockLength;
        this.magnification = magnification;

        if (chosenDevice == null) {
            IJ.log("Looking for the fastest device...");
            System.out.println("Using the fastest device...");
            context = CLContext.create(clPlatformMaxFlop);
            chosenDevice = context.getMaxFlopsDevice();
            IJ.log("Using "+chosenDevice.getName());
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

        clBufferPx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_ONLY);
        clBufferShiftXY = context.createFloatBuffer(2 * nFrameForSRRF, READ_ONLY);
        clBufferGx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gx
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gy
        clBufferGxInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gx
        clBufferGyInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gy
        clBufferOut = context.createFloatBuffer((nReconstructions + 1) * widthM * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);
        clBufferMPmap = context.createFloatBuffer(2 * magnification * magnification, READ_WRITE);

        // Current frame is a 2 element Int buffer:
        // nCurrentFrame[0] is the global current frame in the current SRRF frame (reset every SRRF frame)
        // nCurrentFrame[1] is the local current frame in the current GPU-loaded dataset (reset every tun of the method calculateSRRF (within the gradient calculation))


        // Create the program
        float sigma = fwhm / 2.354f;
        float radius = ((float) ((int) (GxGyMagnification * 2 * sigma))) / GxGyMagnification + 1;    // this reduces the radius for speed, works when using dGauss^4 and 2p+I

        String programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");
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



        int argn;
        argn = 0;
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelInterpolateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGyInt); // make sure type is the same !!

        argn = 0;
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftXY); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!

        kernelIncrementFramePosition.setArg(0, clBufferCurrentFrame); // make sure type is the same !!

        argn = 0;
        kernelCalculateMPmap.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateMPmap.setArg(argn++, clBufferMPmap); // make sure type is the same !!

        argn = 0;
        kernelCalculateStd.setArg(argn++, clBufferOut);


        queue = chosenDevice.createCommandQueue();

        System.out.println("used device memory: " + (  // TODO: add the buffer for MPmap
                        clBufferPx.getCLSize() +
                        clBufferShiftXY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferGxInt.getCLSize() +
                        clBufferGyInt.getCLSize() +
                        clBufferOut.getCLSize() +
                        clBufferCurrentFrame.getCLSize())
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
    public synchronized boolean calculateSRRF(ImageStack imsRawData) { // returns boolean desribing whether it was cancelled by user or not

        assert (imsRawData.getWidth() == width && imsRawData.getHeight() == height);
        int nFrameToLoad = imsRawData.getSize();

//        IJ.log("Uploading raw data to GPU...");
        int id = prof.startTimer();
        fillBuffer(clBufferPx, imsRawData);
        queue.putWriteBuffer(clBufferPx, false);
        prof.recordTime("Uploading data to GPU", prof.endTimer(id));

        // Make kernelCalculateGradient assignment (this kernel also resets the local GPU load frame counter)
//        IJ.log("Calculating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put3DRangeKernel(kernelCalculateGradient, 0, 0, 0, width, height, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

//        IJ.log("Interpolating gradient...");
        id = prof.startTimer();
        queue.finish(); // Make sure everything is done
        queue.put3DRangeKernel(kernelInterpolateGradient, 0, 0, 0, GxGyMagnification * width, GxGyMagnification * height, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelInterpolateGradient", prof.endTimer(id));

        // Make kernelCalculateSRRF assignment
//        IJ.log("Calculating SRRF...");

        int workSize;

        for (int f = 0; f < nFrameToLoad; f++) {

            int nBlocks = widthM * heightM / blockLength + ((widthM * heightM % blockLength == 0) ? 0 : 1);

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

        // Calculate the STD on the OutputArray on the GPU
        queue.put1DRangeKernel(kernelCalculateStd, 0, heightM*heightM,0);

        return false;

    }

    // --- Read the output buffer ---
    public void readSRRFbuffer() {

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferOut, true);
        FloatBuffer bufferSRRF = clBufferOut.getBuffer();

//        ImageStack imsSRRF = new ImageStack(widthM, heightM);
        float[] dataSRRF;
        imsSRRF = new ImageStack(widthM, heightM);

        // Load average
        dataSRRF = new float[widthM * heightM];
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

//        return imsSRRF;
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

    // --- Calculate and get the Macro-pixel map ---
    public ImageStack getMPmap() {

        int id = prof.startTimer();
        queue.put1DRangeKernel(kernelCalculateMPmap,0,2 * magnification * magnification,0);
        prof.recordTime("Calculate MP map", prof.endTimer(id));

        ImageStack imsMPmap = new ImageStack(magnification, 2 * magnification);
        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferMPmap, true);
        FloatBuffer bufferMPmap = clBufferMPmap.getBuffer();

        float[] dataMPmap = new float[2*magnification*magnification];
        for (int n = 0; n < 2*magnification*magnification; n++) {
            if (Float.isNaN(dataMPmap[n])) dataMPmap[n] = 0; // make sure we dont get any weirdness
            dataMPmap[n] = bufferMPmap.get(n);
        }

        imsMPmap.addSlice(new FloatProcessor(magnification, 2 * magnification, dataMPmap));

        return imsMPmap;

    }


}
