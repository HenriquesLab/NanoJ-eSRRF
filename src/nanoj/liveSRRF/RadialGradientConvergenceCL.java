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
    static CLKernel kernelCalculateRGC, kernelCalculateGradient;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, magnification;
    private final float fwhm;

    public int nVectors = 12;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferRGC;

    public RadialGradientConvergenceCL(int width, int height, int magnification, float fwhm) {
        this.width = width;
        this.height = height;
        this.widthM = width * magnification;
        this.heightM = height * magnification;
        this.magnification = magnification;
        this.fwhm = fwhm;

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
        kernelCalculateGradient = program.createCLKernel("calculateGradientRobX");
//        kernelCalculateGradient = program.createCLKernel("calculateGradient");
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
                / 1000000d + "MB");
    }

    public synchronized FloatProcessor calculateRGC(ImageProcessor ip, float shiftX, float shiftY) {
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

        // make kernelCalculateRGC assignment
        argn = 0;
        kernelCalculateRGC.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferRGC); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, magnification ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, fwhm ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, shiftX ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, shiftY ); // make sure type is the same !!

        // asynchronous write of data to GPU device,
        // followed by blocking read to get the computed results back.
        int id = prof.startTimer();
        queue.putWriteBuffer( clBufferPx, false );
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);
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
