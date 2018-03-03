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

    private CLBuffer<FloatBuffer>
            clBufferPx, clBufferPxM,
            clBufferGx, clBufferGy,
            clBufferWM0, clBufferWM1, // weight masks
            clBufferRGC;

    public RadialGradientConvergenceCL(FloatProcessor fpWeightMask, int magnification, float fwhm) {
        this.width = fpWeightMask.getWidth();
        this.height = fpWeightMask.getHeight();
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
        kernelCalculateGradient = program.createCLKernel("calculateGradient");
        kernelCalculateRGC = program.createCLKernel("calculateRadialGradientConvergence");

        clBufferPx = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferGy = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferWM0 = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferWM1 = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferRGC = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);
        clBufferPxM = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);

        // prepare CL buffers
        fillBuffer(clBufferWM0, fpWeightMask);
        queue.putWriteBuffer( clBufferWM0, false ); // copy already the first weight mask into GPU memory

        System.out.println("used device memory: " + (
                        clBufferPx.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferWM0.getCLSize() +
                        clBufferWM1.getCLSize() +
                        clBufferRGC.getCLSize() +
                        clBufferPxM.getCLSize())
                / 1000000d + "MB");
    }

    public synchronized FloatProcessor[] calculateRGC(ImageProcessor ip, float shiftX, float shiftY) {
        assert (ip.getWidth() == width && ip.getHeight() == height);

        FloatProcessor fpPx = ip.convertToFloatProcessor();

        FloatProcessor fpRGC = new FloatProcessor(widthM, heightM);
        FloatProcessor fpPxM = new FloatProcessor(widthM, heightM);

        // prepare CL buffers
        fillBuffer(clBufferPx, fpPx);

        int argn;
        // make kernelCalculateRGC assignment
        argn = 0;
        kernelCalculateGradient.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferWM0 ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferWM1 ); // make sure type is the same !!

        // make kernelCalculateRGC assignment
        argn = 0;
        kernelCalculateRGC.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferWM1 ); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferRGC); // make sure type is the same !!
        kernelCalculateRGC.setArg( argn++, clBufferPxM); // make sure type is the same !!
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
        queue.putReadBuffer(clBufferPxM, true);
        prof.recordTime("RadialGradientConvergence.cl", prof.endTimer(id));

        //clBufferRGC.getBuffer().get((float[]) fpRadiality.getPixels());
        grabBuffer(clBufferRGC, fpRGC);
        grabBuffer(clBufferPxM, fpPxM);

        return new FloatProcessor[] {fpPxM, fpRGC};
    }

    public void release() {
        context.release();
    }

}
