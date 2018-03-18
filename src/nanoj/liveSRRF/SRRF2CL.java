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
import static java.lang.Math.max;
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.grabBuffer;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackFromFloatArray;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;

public class SRRF2CL {

    public boolean DEBUG = false;
    private NanoJProfiler prof = new NanoJProfiler();

    public static String[] reconstructionLabel = new String[]{
            "Raw Average",
            "SRRF Maximum",
            "SRRF Average",
            "SRRF StdDev",
            "SRRF 2nd-Order",
            "SRRF 3rd-Order",
            "SRRF 4th-Order"
    };

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelCalculateRGC, kernelCalculateGradient, kernelCalculateSRRF;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, nFrames, magnification, whM;
    private final float fwhm;



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
        this.whM = widthM * heightM;
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
        clBufferSRRF = context.createFloatBuffer(widthM * heightM * 7, WRITE_ONLY);

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

        FloatBuffer buffer = clBufferSRRF.getBuffer();
        ImageStack imsSRRF = new ImageStack(widthM, heightM);

        // grab interpolated data frame
        float[] SRRF_RAW = new float[whM];
        float SRRF_RAW_MAX = -Float.MAX_VALUE;
        float SRRF_RAW_MIN = Float.MAX_VALUE;
        for(int n=0; n<whM; n++) {
            SRRF_RAW[n] = buffer.get(n);
            SRRF_RAW_MAX = max(SRRF_RAW[n], SRRF_RAW_MAX);
            SRRF_RAW_MIN = min(SRRF_RAW[n], SRRF_RAW_MIN);
        }
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, SRRF_RAW));

        // grab other data frames
        for (int r=1; r<7; r++) {
            int offset = whM * r;

            float[] data = new float[whM];
            float dataMax = -Float.MAX_VALUE;
            float dataMin = Float.MAX_VALUE;

            for(int n=0; n<whM; n++) {
                data[n] = buffer.get(n+offset);
                dataMax = max(data[n], dataMax);
                dataMin = min(data[n], dataMin);
            }

            // renormalise data
            float gain = (SRRF_RAW_MAX - SRRF_RAW_MIN) / dataMax;
            for(int n=0; n<whM; n++) {
                data[n] = (data[n] - dataMin) * gain + SRRF_RAW_MIN;
            }
            imsSRRF.addSlice(new FloatProcessor(widthM, heightM, data));
        }

        return imsSRRF;
    }

    public void release() {
        context.release();
    }
}
