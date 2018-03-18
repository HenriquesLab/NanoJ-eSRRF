package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJProfiler;
import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.grabBuffer;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackFromFloatArray;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;

public class SRRF2CL {

    public boolean DEBUG = false;
    public NanoJProfiler prof = new NanoJProfiler();

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelCalculateRGC, kernelCalculateGradient, kernelCalculateSRRF;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, nFrames, magnification;
    private final float fwhm;

    public int nVectors = 12;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferShiftX, clBufferShiftY,
            clBufferSRRF;

    public SRRF2CL(int width, int height, int nFrames, int magnification, float fwhm) {
        this.nFrames = nFrames;
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
            InputStream programStream = RadialGradientConvergenceCL.class.getResourceAsStream("/SRRF2.cl");
            String programString = IOUtils.toString(programStream);
            programString = programString.replace("MAX_FRAMES", ""+nFrames);
            program = context.createProgram(programString).build();
        } catch (IOException e) {
            e.printStackTrace();
        }
        kernelCalculateGradient = program.createCLKernel("calculateGradientRobX");
//        kernelCalculateGradient = program.createCLKernel("calculateGradient");
        kernelCalculateSRRF = program.createCLKernel("calculateSRRF");

        clBufferPx = context.createFloatBuffer(width * height * nFrames, READ_ONLY);
        clBufferShiftX = context.createFloatBuffer(nFrames, READ_ONLY);
        clBufferShiftY = context.createFloatBuffer(nFrames, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height * nFrames, READ_WRITE);
        clBufferGy = context.createFloatBuffer(width * height * nFrames, READ_WRITE);
        clBufferSRRF = context.createFloatBuffer(widthM * heightM * 10, WRITE_ONLY);

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        //clBufferRGC.getCLSize() +
                        clBufferSRRF.getCLSize())
                / 1000000d + "MB");
    }

    public synchronized ImageStack calculateSRRF(ImageStack ims, float[] shiftX, float[] shiftY) {
        assert (ims.getWidth() == width && ims.getHeight() == height);
        assert (ims.getSize() <= nFrames);

        float[] dataPx = ImageStackToFloatArray(ims);
        int nFrames = ims.getSize();
        int argn;

        // prepare and upload CL buffers
        IJ.showStatus("Uploading data to GPU...");
        int id = prof.startTimer();
        fillBuffer(clBufferPx, dataPx);
        fillBuffer(clBufferShiftX, shiftX);
        fillBuffer(clBufferShiftY, shiftY);
        queue.putWriteBuffer( clBufferPx, false );
        queue.putWriteBuffer( clBufferShiftX, false );
        queue.putWriteBuffer( clBufferShiftY, false );
        prof.recordTime("uploading data to GPU", prof.endTimer(id));

        // make kernelCalculateRGC assignment
        argn = 0;
        IJ.showStatus("Calculating gradient...");
        id = prof.startTimer();
        kernelCalculateGradient.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, nFrames ); // make sure type is the same !!
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

        // make kernelCalculateSRRF assignment
        argn = 0;
        IJ.showStatus("Calculating SRRF...");
        id = prof.startTimer();
        kernelCalculateSRRF.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferShiftX ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferShiftY ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferSRRF); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, nFrames); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, magnification ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, fwhm ); // make sure type is the same !!
        queue.put2DRangeKernel(kernelCalculateSRRF, 0, 0, widthM, heightM, 0, 0);
        IJ.showStatus("Waiting for GPU to finish...");
        queue.finish();
        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

        // download CL buffers
        IJ.showStatus("Downloading data from GPU...");
        id = prof.startTimer();
        queue.putReadBuffer(clBufferSRRF, true);
        float[] data = new float[widthM * heightM * 10];
        grabBuffer(clBufferSRRF, data);
        prof.recordTime("dowloading data from GPU", prof.endTimer(id));

        ImageStack imsSRRF = ImageStackFromFloatArray(data, widthM, heightM);
        return imsSRRF;
    }

    public void release() {
        context.release();
    }
}
