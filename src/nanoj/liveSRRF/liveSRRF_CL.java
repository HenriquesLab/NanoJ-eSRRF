package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.getResourceAsString;
import static nanoj.core2.NanoJCL.replaceFirst;

public class liveSRRF_CL {

    // Basic formats
    private int width,
            height,
            widthM,
            heightM;

    private final int GxGyMagnification = 2;
    private final float vxy_offset = 0.5f;
    private final int vxy_ArrayShift = 1;

    private final int nReconstructions = 2; // Currently only STD and AVG

    private final String[] reconstructionLabels = {"AVG", "STD", "RawInt"};

    // Advanced formats
    private NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static CLContext context;
    static CLProgram programLiveSRRF;
    static CLKernel kernelCalculateGradient,
            kernelInterpolateGradient,
            kernelCalculateSRRF,
            kernelIncrementFramePosition,
            kernelResetFramePosition;

    static CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferGxInt, clBufferGyInt,
            clBufferShiftX, clBufferShiftY,
            clBufferOut;

    private CLBuffer<IntBuffer>
            clBufferCurrentFrame;

    // --- Initialization method ---
    public liveSRRF_CL(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF) {

        this.width = width;
        this.height = height;
        this.heightM = height * magnification;
        this.widthM = width * magnification;

        context = CLContext.create();
        System.out.println("created " + context);

        // Select fastest device
        CLDevice device = context.getMaxFlopsDevice();
        System.out.println("using " + device);

        clBufferPx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_ONLY);
        clBufferShiftX = context.createFloatBuffer(nFrameForSRRF, READ_ONLY);
        clBufferShiftY = context.createFloatBuffer(nFrameForSRRF, READ_ONLY);
        clBufferGx = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gx
        clBufferGy = context.createFloatBuffer(nFramesOnGPU * width * height, READ_WRITE); // single frame Gy
        clBufferGxInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gx
        clBufferGyInt = context.createFloatBuffer(nFramesOnGPU * 4 * width * height, READ_WRITE); // single frame Gy
        clBufferOut = context.createFloatBuffer((nReconstructions + 1) * widthM * nReconstructions * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(1, READ_WRITE);

        // Create the program
        float sigma = fwhm / 2.354f;
        String programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");
        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
        programString = replaceFirst(programString, "$FWHM$", "" + fwhm);
        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);
        programString = replaceFirst(programString, "$GXGYMAGNIFICATION$", "" + GxGyMagnification);

        programString = replaceFirst(programString, "$SIGMA$", "" + sigma);
        programString = replaceFirst(programString, "$RADIUS$", "" + ((int) (sigma * 2) + 1));
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

        programString = replaceFirst(programString, "$NFRAMEONGPU$", "" + nFramesOnGPU);


        programLiveSRRF = context.createProgram(programString).build();

        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");
        kernelInterpolateGradient = programLiveSRRF.createCLKernel("calculateGradientInterpolation");
        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateRadialGradientConvergence");

        kernelIncrementFramePosition = programLiveSRRF.createCLKernel("kernelIncrementFramePosition");
        kernelResetFramePosition = programLiveSRRF.createCLKernel("kernelResetFramePosition");

        queue = device.createCommandQueue();

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferShiftX.getCLSize() +
                        clBufferShiftY.getCLSize() +
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
        int id = prof.startTimer();
        fillBuffer(clBufferShiftX, shiftX); // TODO: load ShiftXY as a single buffer??
        queue.putWriteBuffer(clBufferShiftX, false);
        fillBuffer(clBufferShiftY, shiftY);
        queue.putWriteBuffer(clBufferShiftY, false);
        prof.recordTime("Uploading shift arrays to GPU", prof.endTimer(id));
        queue.finish(); // Make sure everything is done... // TODO: do we need to do this for every method?

    }

    // --- Release GPU context ---
    public void release() {
        context.release();
    }

    // --- Reset the SRRF frame counter ---
    public void resetFramePosition() {
        int id = prof.startTimer();
        kernelResetFramePosition.setArg(0, clBufferCurrentFrame); // make sure type is the same !!
        queue.put1DRangeKernel(kernelResetFramePosition, 0, 1, 1);
        prof.recordTime("Reset SRRF frame counter", prof.endTimer(id));
    }


    // --- Calculate SRRF images ---
    public synchronized void calculateSRRF(ImagePlus imp, int indexStart, int nFrameToLoad) {

        int argn;
        assert (imp.getWidth() == width && imp.getHeight() == height);
        ImageStack imsRawData = new ImageStack(width, height);

        for (int f = 0; f < nFrameToLoad; f++) {
            imp.setSlice(indexStart + f);
            imsRawData.addSlice(imp.getProcessor());
        }

        IJ.log("----------------------------");
        IJ.log("Uploading raw data to GPU...");
        int id = prof.startTimer();
        fillBuffer(clBufferPx, imsRawData);
        queue.putWriteBuffer(clBufferPx, false);
        prof.recordTime("Uploading data to GPU", prof.endTimer(id));

        // Make kernelCalculateGradient assignment
        IJ.log("Calculating gradient...");
        id = prof.startTimer();
        argn = 0;
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        queue.put3DRangeKernel(kernelCalculateGradient, 0, 0, 0, width, height, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

        IJ.log("Interpolating gradient...");
        id = prof.startTimer();
        argn = 0;
        kernelInterpolateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelInterpolateGradient.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        queue.put3DRangeKernel(kernelInterpolateGradient, 0, 0, 0, GxGyMagnification * width, GxGyMagnification * height, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelInterpolateGradient", prof.endTimer(id));

        // Make kernelCalculateSRRF assignment
        argn = 0;
        IJ.log("Calculating SRRF...");
        id = prof.startTimer();
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGxInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGyInt); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferOut); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftX); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftY); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferCurrentFrame); // make sure type is the same !!
        queue.put3DRangeKernel(kernelCalculateSRRF, 0, 0, 0, widthM, heightM, nFrameToLoad, 0, 0, 0);
        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

        id = prof.startTimer();
        kernelIncrementFramePosition.setArg(0, clBufferCurrentFrame); // make sure type is the same !!
        queue.put1DRangeKernel(kernelIncrementFramePosition, 0, 2, 2);
        prof.recordTime("Increment frame count", prof.endTimer(id));

        queue.finish(); // Make sure everything is done... //TODO: do we need to do this for every method?
    }

    public ImageStack readSRRFbuffer() {

        queue.putReadBuffer(clBufferOut, true);
        FloatBuffer bufferSRRF = clBufferOut.getBuffer();

        ImageStack imsSRRF = new ImageStack(widthM, heightM);
        float[] dataSRRF = new float[widthM * heightM];

        for (int r = 0; r < (nReconstructions + 1); r++) {
            for (int n = 0; n < widthM * heightM; n++) {
                dataSRRF[n] = bufferSRRF.get(widthM * heightM + n);
                if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we dont get any weirdness
            }

            imsSRRF.addSlice(reconstructionLabels[r], new FloatProcessor(widthM, heightM, dataSRRF));
        }

        return imsSRRF;
    }


}
