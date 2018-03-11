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

public class CompareGradientEstimationCL {

    public NanoJProfiler prof = new NanoJProfiler();

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelCalculateGradient, kernelInterpolateGradient;
    static CLCommandQueue queue;

    private final int width, height;
    private float vxy_offset;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferGxInt, clBufferGyInt;


    public CompareGradientEstimationCL(int width, int height, String GradMethod) {
        this.width = width;
        this.height = height;

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
            program = context.createProgram(CompareGradientEstimationCL.class.getResourceAsStream("/RadialGradientConvergence.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }

        // Choose gradient kernel
        if      (GradMethod.equals("3-point gradient (classic)")){
            kernelCalculateGradient = program.createCLKernel("calculateGradient");
            // in this method the gradient is calculated in the centre of the designated pixel
            this.vxy_offset = 0.5f;}

        else if (GradMethod.equals("Robert's cross local gradient")){
            kernelCalculateGradient = program.createCLKernel("calculateGradientRobX");
            // in this method the gradient is estimated at the crossing of the pixels, therefore at an offset of 0.5 with respect to the pixel centre
            this.vxy_offset = 0.0f;}

        else if (GradMethod.equals("2-point local + interpolation")){
            kernelCalculateGradient = program.createCLKernel("calculateGradient_2point");
            kernelInterpolateGradient = program.createCLKernel("calculateGradient2p_Interpolation");

            // this method currently wrecks havock on your data
            this.vxy_offset = 0.0f;
            clBufferGxInt = context.createFloatBuffer(4*width * height, READ_WRITE);
            clBufferGyInt = context.createFloatBuffer(4*width * height, READ_WRITE);}


        clBufferPx = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height, READ_WRITE);
        clBufferGy = context.createFloatBuffer(width * height, READ_WRITE);

        System.out.println("used device memory: " + (
                        clBufferPx.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize())
                / 1000000d + "MB"); // TODO: add size of GxInt and GyInt in total buffer size
    }

    public synchronized FloatProcessor[] calculateGxGy(ImageProcessor ip, float shiftX, float shiftY, String GradChosenMethod) {

        int GxGyWidth;
        int GxGyHeight;

//        if (GradChosenMethod.equals("2-point local + interpolation")){
//            GxGyWidth = 2*width;
//            GxGyHeight = 2*height;}
//        else {

            GxGyWidth = width;
            GxGyHeight = height;
//    }

        assert (ip.getWidth() == GxGyWidth && ip.getHeight() == GxGyHeight);

        FloatProcessor fpPx = ip.convertToFloatProcessor();
        FloatProcessor fpGx = new FloatProcessor(GxGyWidth, GxGyHeight);
        FloatProcessor fpGy = new FloatProcessor(GxGyWidth, GxGyHeight);

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



        // asynchronous write of data to GPU device,
        // followed by blocking read to get the computed results back.
        int id = prof.startTimer();
        queue.putWriteBuffer( clBufferPx, false );
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);

        if (GradChosenMethod.equals("2-point local + interpolation")){
            queue.put2DRangeKernel(kernelInterpolateGradient, 0, 0, 2*width, 2*height, 0, 0);
            queue.finish();
            queue.putReadBuffer(clBufferGxInt, true);
            queue.putReadBuffer(clBufferGyInt, true);

            prof.recordTime("RadialGradientConvergence.cl", prof.endTimer(id));
            grabBuffer(clBufferGxInt, fpGx);
            grabBuffer(clBufferGyInt, fpGy);
        }
        else {
            queue.finish();
            queue.putReadBuffer(clBufferGx, true);
            queue.putReadBuffer(clBufferGy, true);
            prof.recordTime("RadialGradientConvergence.cl", prof.endTimer(id));

            grabBuffer(clBufferGx, fpGx);
            grabBuffer(clBufferGy, fpGy);
        }


        FloatProcessor[] fpGxy = new FloatProcessor[2];
        fpGxy[0] = fpGx;
        fpGxy[1] = fpGy;

        return fpGxy;
    }

    public void release() {
        context.release();
    }

}
