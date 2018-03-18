package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;
import org.apache.commons.io.IOUtils;

import java.io.IOException;
import java.io.InputStream;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.grabBuffer;
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
            "SRRF 4th-Order",
            "SRRF Fusion"
    };

    static CLContext context;
    static CLProgram programSRRF2, programConvolve2DIntegratedGaussian;
    static CLKernel kernelCalculateGradient, kernelCalculateSRRF, kernelConvolveH, kernelConvolveV;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, nFrames, magnification, whM;
    private final float fwhm;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferShiftX, clBufferShiftY,
            clBufferSRRF, clBufferSRRF_CVH, clBufferSRRF_CV;

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
            programSRRF2 = context.createProgram(programString).build();

            programConvolve2DIntegratedGaussian = context.createProgram(RadialGradientConvergenceCL.class.getResourceAsStream("/Convolve2DIntegratedGaussian.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }
        kernelCalculateGradient = programSRRF2.createCLKernel("calculateGradientRobX");
        kernelCalculateSRRF = programSRRF2.createCLKernel("calculateSRRF");
        kernelConvolveH = programConvolve2DIntegratedGaussian.createCLKernel("convolveHorizontal");
        kernelConvolveV = programConvolve2DIntegratedGaussian.createCLKernel("convolveVertical");

        clBufferPx = context.createFloatBuffer(width * height * nFrames, READ_ONLY);
        clBufferShiftX = context.createFloatBuffer(nFrames, READ_ONLY);
        clBufferShiftY = context.createFloatBuffer(nFrames, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height * nFrames, READ_WRITE);
        clBufferGy = context.createFloatBuffer(width * height * nFrames, READ_WRITE);
        clBufferSRRF = context.createFloatBuffer(widthM * heightM * 7, READ_WRITE);
        clBufferSRRF_CVH = context.createFloatBuffer(widthM * heightM * 7, READ_WRITE);
        clBufferSRRF_CV  = context.createFloatBuffer(widthM * heightM * 7, READ_WRITE);

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferShiftX.getCLSize() +
                        clBufferShiftY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferSRRF.getCLSize() +
                        clBufferSRRF_CVH.getCLSize() +
                        clBufferSRRF_CV.getCLSize())
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
        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

        argn = 0;
        IJ.showStatus("Convolving SRRF...");
        float sigmaM = magnification * fwhm / 2.354f;
        id = prof.startTimer();
        kernelConvolveH.setArg( argn++, clBufferSRRF ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, clBufferSRRF_CVH ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, 0 ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, 7 ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, sigmaM); // make sure type is the same !!
        queue.put2DRangeKernel(kernelConvolveH, 0, 0, widthM, heightM, 0, 0);
        argn = 0;
        kernelConvolveV.setArg( argn++, clBufferSRRF_CVH ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, clBufferSRRF_CV ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, 0 ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, 7 ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, sigmaM ); // make sure type is the same !!
        queue.put2DRangeKernel(kernelConvolveV, 0, 0, widthM, heightM, 0, 0);
        prof.recordTime("kernelConvolveH & kernelConvolveV", prof.endTimer(id));

        id = prof.startTimer();
        queue.finish(); // make sure everything is done...
        prof.recordTime("waiting for queue to finish", prof.endTimer(id));

        // download CL buffers
        IJ.showStatus("Downloading data from GPU...");
        id = prof.startTimer();
        queue.putReadBuffer(clBufferSRRF, true);
        queue.putReadBuffer(clBufferSRRF_CV, true);
        prof.recordTime("downloading data from GPU", prof.endTimer(id));

        FloatBuffer bufferSRRF = clBufferSRRF.getBuffer();
        FloatBuffer bufferSRRF_CV = clBufferSRRF_CV.getBuffer();
        ImageStack imsSRRF = new ImageStack(widthM, heightM);

        // grab interpolated data frame and calculate its max and min
        IJ.showStatus("Preparing data...");
        float[] RAW_AVE = new float[whM];
        float RAW_AVE_MAX = -Float.MAX_VALUE;
        float RAW_AVE_MIN =  Float.MAX_VALUE;
        for(int n=0; n<whM; n++) {
            RAW_AVE[n] = bufferSRRF.get(n);
            RAW_AVE_MAX = max(RAW_AVE[n], RAW_AVE_MAX);
            RAW_AVE_MIN = min(RAW_AVE[n], RAW_AVE_MIN);
        }

        // grab data from convolved frames
        double[][] errorMap = new double[7][whM];
        double[] errorMax = new double[whM];

        for (int r=0; r<7; r++) {
            int offset = whM * r;

            float[] data = new float[whM];
            float[] dataConvolved = new float[whM];
            float dataMax = -Float.MAX_VALUE;
            float dataMin = Float.MAX_VALUE;

            for(int n=0; n<whM; n++) {
                data[n] = bufferSRRF.get(n+offset);
                dataConvolved[n] = bufferSRRF_CV.get(n+offset);
                dataMax = max(dataConvolved[n], dataMax);
                dataMin = min(dataConvolved[n], dataMin);
            }

            // renormalise data and calculate erro map
            float gain = (RAW_AVE_MAX - RAW_AVE_MIN) / dataMax;
            if (r>0) for (int n = 0; n < whM; n++) data[n] = (data[n] - dataMin) * gain + RAW_AVE_MIN;

            for (int n = 0; n < whM; n++) {
                dataConvolved[n] = (dataConvolved[n] - dataMin) * gain + RAW_AVE_MIN;
                errorMap[r][n] = abs(RAW_AVE[n] - dataConvolved[n]);
                errorMax[n] = max(errorMax[n], errorMap[r][n]);
            }

            imsSRRF.addSlice(new FloatProcessor(widthM, heightM, data));
        }

        IJ.showStatus("Calculating SRRF Fusion...");
        float[] fusion = new float[whM];
        double[] weightSum = new double[whM];

        for (int r=0; r<7; r++) {
            float[] pixelsSRRF = (float[]) imsSRRF.getProcessor(r+1).getPixels();
            for (int n = 0; n < whM; n++) {
                double w = (errorMax[n]) / max(errorMap[r][n], 1);
                fusion[n] += pixelsSRRF[n] * w;
                weightSum[n] += w;
            }
        }
        for (int n = 0; n < whM; n++) fusion[n] /= weightSum[n];
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, fusion));

        return imsSRRF;
    }

    public void release() {
        context.release();
    }
}
