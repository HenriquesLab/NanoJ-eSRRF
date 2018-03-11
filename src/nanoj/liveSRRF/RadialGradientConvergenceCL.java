package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJProfiler;

import java.io.IOException;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.System.nanoTime;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.grabBuffer;

public class RadialGradientConvergenceCL {

    public boolean DEBUG = false;
    public NanoJProfiler prof = new NanoJProfiler();

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelCalculateRGC, kernelCalculateGradient, kernelInterpolateGradient;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, magnification;
    private final float fwhm;
    private float vxy_offset;
//    private String GradMethod;

    public int nVectors = 12;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferGxInt, clBufferGyInt,
            clBufferRGC;


    public RadialGradientConvergenceCL(int width, int height, int magnification, float fwhm, String GradMethod) {
        this.width = width;
        this.height = height;
        this.widthM = width * magnification;
        this.heightM = height * magnification;
        this.magnification = magnification;
        this.fwhm = fwhm;
//        this.GradMethod = GradMethod;

        //Create a context from GPU
        //context = CLContext.create( CLDevice.Type.GPU );
        context = CLContext.create();
        System.out.println("created " + context);

        // Select fastest device
        CLDevice device = context.getMaxFlopsDevice();
        System.out.println("using " + device);

        // create command queue on device.
        queue = device.createCommandQueue();

        // create the program
        try {
            program = context.createProgram(RadialGradientConvergenceCL.class.getResourceAsStream("/RadialGradientConvergence.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Choose gradient kernel
        if      (GradMethod.equals("3-point gradient (classic)")){
            kernelCalculateGradient = program.createCLKernel("calculateGradient");
            // in this method the gradient is calculated in the centre of the designated pixel
            this.vxy_offset = 0.5f; // signed distance between the (0,0) continuous space position and the position at which the vector in the (0,0) position in the Gradient image is calculated
        }

        else if (GradMethod.equals("Robert's cross local gradient")){
            kernelCalculateGradient = program.createCLKernel("calculateGradientRobX");
            // in this method the gradient is estimated at the crossing of the pixels, therefore at an offset of 0.5 with respect to the pixel centre
            this.vxy_offset = 0.0f;
        }

        else if (GradMethod.equals("2-point local + interpolation")){
            kernelCalculateGradient = program.createCLKernel("calculateGradient_2point");
            kernelInterpolateGradient = program.createCLKernel("calculateGradient2p_Interpolation");

            // this method currently wrecks havock on your data
            this.vxy_offset = 0.0f;
            clBufferGxInt = context.createFloatBuffer(4*width * height, READ_WRITE);
            clBufferGyInt = context.createFloatBuffer(4*width * height, READ_WRITE);
        }


        // Prepare kernel for GRC calculation
        kernelCalculateRGC = program.createCLKernel("calculateRadialGradientConvergence");

        clBufferPx = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferGy = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferRGC = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);

        System.out.println("used device memory: " + (
                        clBufferPx.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferRGC.getCLSize())
                / 1000000d + "MB"); // TODO: add the clBufferGxInt and GyInt to the buffer size calculation
    }

    public synchronized FloatProcessor calculateRGC(ImageProcessor ip, float shiftX, float shiftY, String GradChosenMethod) {
        assert (ip.getWidth() == width && ip.getHeight() == height);

        FloatProcessor fpPx = ip.convertToFloatProcessor();

        FloatProcessor fpRGC = new FloatProcessor(widthM, heightM);
        FloatProcessor fpInterpolated = new FloatProcessor(widthM, heightM);

        // prepare CL buffers
        fillBuffer(clBufferPx, fpPx);

        int argn;

        // make kernelCalculateRGC assignment
        argn = 0;
        kernelCalculateGradient.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGy ); // make sure type is the same !!

        argn = 0;
        if (GradChosenMethod.equals("2-point local + interpolation")){
            kernelInterpolateGradient.setArg( argn++, clBufferGx ); // make sure type is the same !!
            kernelInterpolateGradient.setArg( argn++, clBufferGy ); // make sure type is the same !!
            kernelInterpolateGradient.setArg( argn++, clBufferGxInt ); // make sure type is the same !!
            kernelInterpolateGradient.setArg( argn++, clBufferGyInt ); // make sure type is the same !!
        }

        // make kernelCalculateRGC assignment
        argn = 0;
        kernelCalculateRGC.setArg( argn++, clBufferPx ); // make sure type is the same !!
        if (GradChosenMethod.equals("2-point local + interpolation")){
            kernelCalculateRGC.setArg( argn++, clBufferGxInt ); // make sure type is the same !!
            kernelCalculateRGC.setArg( argn++, clBufferGyInt ); // make sure type is the same !!
            }
        else {
            kernelCalculateRGC.setArg( argn++, clBufferGx ); // make sure type is the same !!
            kernelCalculateRGC.setArg( argn++, clBufferGy ); // make sure type is the same !!
            }

        kernelCalculateRGC.setArg( argn++, clBufferRGC); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, magnification ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, fwhm ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, shiftX ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, shiftY ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, vxy_offset ); // make sure type is the same !!

        // asynchronous write of data to GPU device,
        // followed by blocking read to get the computed results back.
        int id = prof.startTimer();
        queue.putWriteBuffer( clBufferPx, false );
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);

        if (GradChosenMethod.equals("2-point local + interpolation")){
            queue.put2DRangeKernel(kernelInterpolateGradient, 0, 0, 2*width, 2*height, 0, 0);
        }

        queue.put2DRangeKernel(kernelCalculateRGC, 0, 0, widthM, heightM, 0, 0);
        queue.finish();
        queue.putReadBuffer(clBufferRGC, true);
        prof.recordTime("RadialGradientConvergence.cl", prof.endTimer(id));

        grabBuffer(clBufferRGC, fpRGC);

        return fpRGC;
    }

    public void release() {
        context.release();
    }

}
