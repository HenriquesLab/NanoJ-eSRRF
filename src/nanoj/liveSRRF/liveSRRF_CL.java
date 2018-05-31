package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import nanoj.core2.NanoJProfiler;
import scala.Int;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static nanoj.core2.NanoJCL.getResourceAsString;
import static nanoj.core2.NanoJCL.replaceFirst;

public class liveSRRF_CL {

    //Basic formats
    private int width, height, magnification, sensitivity, widthM, heightM;
    private float fwhm;

    // Advanced formats
    private NanoJProfiler prof = new NanoJProfiler();

    // OpenCL formats
    static CLContext context;
    static CLProgram programLiveSRRF;
    static CLKernel kernelCalculateGradient, kernelInterpolateGradient, kernelCalculateSRRF, kernelIncrementFramePosition;
    static CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferGxInt, clBufferGyInt,
            clBufferShiftX, clBufferShiftY,
            clBufferSRRFavg, clBufferSRRFstd,
            clBufferRawInterpolated;

    private CLBuffer<IntBuffer>
            clBufferCurrentFrame;


    public liveSRRF_CL(int width, int height, int magnification, float fwhm, int sensitivity, int nFramesOnGPU, int nFrameForSRRF) {

        this.width = width;
        this.height = height;
        this.magnification = magnification;
        this.fwhm = fwhm;
        this.sensitivity = sensitivity;
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
        clBufferGx = context.createFloatBuffer(width * height, READ_WRITE); // single frame Gx
        clBufferGy = context.createFloatBuffer(width * height, READ_WRITE); // single frame Gy
        clBufferGxInt = context.createFloatBuffer(4 * width * height, READ_WRITE); // single frame Gx
        clBufferGyInt = context.createFloatBuffer(4 * width * height, READ_WRITE); // single frame Gy
        clBufferSRRFavg = context.createFloatBuffer(widthM * heightM, WRITE_ONLY); // single frame cumulative AVG projection of RGC
        clBufferSRRFstd = context.createFloatBuffer(widthM * heightM, WRITE_ONLY); // single frame cumulative STD projection of RGC
        clBufferRawInterpolated = context.createFloatBuffer(widthM * heightM, WRITE_ONLY); // single frame cumulative STD projection of RGC
        clBufferCurrentFrame = context.createIntBuffer(2, READ_WRITE);

        // Create the program
        float sigma = fwhm / 2.354f;
        String programString = getResourceAsString(liveSRRF_CL.class, "liveSRRF.cl");
        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
        programString = replaceFirst(programString, "$FWHM$", "" + fwhm);
        programString = replaceFirst(programString, "$SENSITIVITY$", "" + sensitivity);

        programString = replaceFirst(programString, "$SIGMA$", "" + sigma);
        programString = replaceFirst(programString, "$RADIUS$", "" + ((int) (sigma * 2) + 1));
        programString = replaceFirst(programString, "$WIDTH$", "" + width);
        programString = replaceFirst(programString, "$HEIGHT$", "" + height);
        programString = replaceFirst(programString, "$WH$", "" + (width * height));
        programString = replaceFirst(programString, "$WM$", "" + (width * magnification));
        programString = replaceFirst(programString, "$HM$", "" + (height * magnification));
        programString = replaceFirst(programString, "$WHM$", "" + (width * height * magnification * magnification));
        programLiveSRRF = context.createProgram(programString).build();

        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");
        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateRadialGradientConvergence");
        kernelInterpolateGradient = programLiveSRRF.createCLKernel("calculateGradientInterpolation");
        kernelIncrementFramePosition = programLiveSRRF.createCLKernel("kernelIncrementFramePosition");


        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferShiftX.getCLSize() +
                        clBufferShiftY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferGxInt.getCLSize() +
                        clBufferGyInt.getCLSize() +
                        clBufferSRRFavg.getCLSize() +
                        clBufferSRRFstd.getCLSize() +
                        clBufferRawInterpolated.getCLSize() +
                        clBufferCurrentFrame.getCLSize())
                / 1000000d + "MB");
    }


}
