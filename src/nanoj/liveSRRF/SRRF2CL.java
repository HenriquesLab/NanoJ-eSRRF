package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;
import nanoj.core2.NanoJThreadExecutor;

import java.io.IOException;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static nanoj.core2.NanoJCL.fillBuffer;
import static nanoj.core2.NanoJCL.getResourceAsString;
import static nanoj.core2.NanoJCL.replaceFirst;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;

public class SRRF2CL {

    public boolean DEBUG = false;
    private NanoJProfiler prof = new NanoJProfiler();

    public static String[] reconstructionLabel = new String[]{
            "Raw Average",
            "SRRF2 Maximum",
            "SRRF2 Average",
            "SRRF2 StdDev",
            "SRRF2 Accum. TL1",
            "SRRF2 Accum. TL2",
            "SRRF2 Accum. TL3",
            "SRRF2 Accum. TL4",
            "SRRF2 Accum. TL5",
            "SRRF2 Fusion"
    };

    static CLContext context;
    static CLProgram programSRRF2, programConvolve2DIntegratedGaussian;
    static CLKernel kernelCalculateGradient, kernelCalculateSRRF, kernelConvolveH, kernelConvolveV;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, nFrames, magnification, whM;
    private final float fwhm;
    private int nGPUReconstructions = reconstructionLabel.length-1;

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
            float sigma = fwhm/2.354f;
            String programString = getResourceAsString(RadialGradientConvergenceCL.class, "SRRF2.cl");
            programString = replaceFirst(programString,"$MAX_FRAMES$", ""+nFrames);
            programString = replaceFirst(programString,"$MAGNIFICATION$", ""+magnification);
            programString = replaceFirst(programString,"$FWHM$", ""+fwhm);
            programString = replaceFirst(programString,"$SIGMA$", ""+sigma);
            programString = replaceFirst(programString,"$WIDTH$", ""+width);
            programString = replaceFirst(programString,"$HEIGHT$", ""+height);
            programString = replaceFirst(programString,"$WH$", ""+(width*height));
            programString = replaceFirst(programString,"$WM$", ""+(width*magnification));
            programString = replaceFirst(programString,"$HM$", ""+(height*magnification));
            programString = replaceFirst(programString,"$WHM$", ""+(width*height*magnification*magnification));
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
        clBufferSRRF = context.createFloatBuffer(widthM * heightM * nGPUReconstructions, READ_WRITE);
        clBufferSRRF_CVH = context.createFloatBuffer(widthM * heightM * nGPUReconstructions, READ_WRITE);
        clBufferSRRF_CV  = context.createFloatBuffer(widthM * heightM * nGPUReconstructions, WRITE_ONLY);

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
        showStatus("Uploading data to GPU...");
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
        showStatus("Calculating gradient...");
        id = prof.startTimer();
        kernelCalculateGradient.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateGradient.setArg( argn++, nFrames ); // make sure type is the same !!
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

        // make kernelCalculateSRRF assignment
        argn = 0;
        showStatus("Calculating SRRF...");
        id = prof.startTimer();
        kernelCalculateSRRF.setArg( argn++, clBufferPx ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferGx ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferGy ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferShiftX ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferShiftY ); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, clBufferSRRF); // make sure type is the same !!
        kernelCalculateSRRF.setArg( argn++, nFrames); // make sure type is the same !!
        queue.put2DRangeKernel(kernelCalculateSRRF, 0, 0, widthM, heightM, 0, 0);
        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));

        argn = 0;
        showStatus("Convolving SRRF...");
        float sigmaM = magnification * fwhm / 2.354f;
        id = prof.startTimer();
        kernelConvolveH.setArg( argn++, clBufferSRRF ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, clBufferSRRF_CVH ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, 0 ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, nGPUReconstructions ); // make sure type is the same !!
        kernelConvolveH.setArg( argn++, sigmaM); // make sure type is the same !!
        queue.put2DRangeKernel(kernelConvolveH, 0, 0, widthM, heightM, 0, 0);
        argn = 0;
        kernelConvolveV.setArg( argn++, clBufferSRRF_CVH ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, clBufferSRRF_CV ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, 0 ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, nGPUReconstructions ); // make sure type is the same !!
        kernelConvolveV.setArg( argn++, sigmaM ); // make sure type is the same !!
        queue.put2DRangeKernel(kernelConvolveV, 0, 0, widthM, heightM, 0, 0);
        prof.recordTime("kernelConvolveH & kernelConvolveV", prof.endTimer(id));

        id = prof.startTimer();
        queue.finish(); // make sure everything is done...
        prof.recordTime("waiting for queue to finish", prof.endTimer(id));

        // download CL buffers
        showStatus("Downloading data from GPU...");
        id = prof.startTimer();
        queue.putReadBuffer(clBufferSRRF, true);
        queue.putReadBuffer(clBufferSRRF_CV, true);
        prof.recordTime("downloading data from GPU", prof.endTimer(id));

        FloatBuffer bufferSRRF = clBufferSRRF.getBuffer();
        FloatBuffer bufferSRRF_CV = clBufferSRRF_CV.getBuffer();
        ImageStack imsSRRF = new ImageStack(widthM, heightM);

        // grab interpolated data frame and calculate its max and min
        showStatus("Preparing data...");
        float[] RAW_AVE = new float[whM];
        for(int n=0; n<whM; n++) RAW_AVE[n] = bufferSRRF.get(n);

        // grab data from convolved frames
        float[][] errorMap = new float[nGPUReconstructions][whM];
        float[] errorMax = new float[whM];
        NanoJThreadExecutor NTE = new NanoJThreadExecutor(true);
        ThreadedReadBufferAndNormalise[] threadList = new ThreadedReadBufferAndNormalise[nGPUReconstructions];

        for (int r=0; r<nGPUReconstructions; r++) {
            int offset = whM * r;
            // this thread will read the buffer for SRRF and SRRF_CV plus rescale their intensity to match RAW_AVE
            ThreadedReadBufferAndNormalise t = new ThreadedReadBufferAndNormalise(RAW_AVE, bufferSRRF, bufferSRRF_CV, offset);
            NTE.execute(t);
            threadList[r] = t;
        }
        NTE.finish();

        // calculate the error map and send off each reconstruction to imsSRRF
        for (int r=0; r<nGPUReconstructions; r++) {
            float[] dataSRRF = threadList[r].dataSRRF;
            errorMap[r] = threadList[r].errorMap;
            errorMax[r] = threadList[r].errorMax;
            imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));
        }

        showStatus("Calculating SRRF Fusion...");
        float[] fusion = new float[whM];
        double[] weightSum = new double[whM];

        for (int r=0; r<7; r++) {
            float[] pixelsSRRF = (float[]) imsSRRF.getProcessor(r+1).getPixels();
            for (int n = 0; n < whM; n++) {
//                double w = (errorMax[n]) / max(errorMap[r][n], 1);
                double w = 1 / max(errorMap[r][n], 1);
                fusion[n] += pixelsSRRF[n] * w;
                weightSum[n] += w;
            }
        }
        for (int n = 0; n < whM; n++) fusion[n] /= weightSum[n];
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, fusion));

        return imsSRRF;
    }

    public void showStatus(String text) {
        IJ.showStatus("SRRF2: "+text);
    }

    public void release() {
        context.release();
    }

    class ThreadedReadBufferAndNormalise extends Thread {
        private final float[] dataRef;
        private final FloatBuffer bufferSRRF;
        private final FloatBuffer bufferSRRF_CV;
        private int offset;
        public float g, o;
        public float[] dataSRRF, dataSRRF_CV, errorMap;
        public float errorMax = - Float.MAX_VALUE;

        public ThreadedReadBufferAndNormalise(float[] dataRef, FloatBuffer bufferSRRF, FloatBuffer bufferSRRF_CV, int offset) {
            this.dataRef = dataRef;
            this.bufferSRRF = bufferSRRF;
            this.bufferSRRF_CV = bufferSRRF_CV;
            this.offset = offset;
        }

        @Override
        public void run() {

            int nPixels = dataRef.length;
            dataSRRF = new float[nPixels];
            dataSRRF_CV = new float[nPixels];
            errorMap = new float[nPixels];

            for (int n=0; n<nPixels; n++) {
                dataSRRF[n] = bufferSRRF.get(n+offset);
                dataSRRF_CV[n] = bufferSRRF_CV.get(n+offset);
                if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we dont get any weirdness
                if (Float.isNaN(dataSRRF_CV[n])) dataSRRF_CV[n] = 0; // make sure we dont get any weirdness
            }

            // now calculate linear regression
            float[] x = dataSRRF_CV;
            float[] y = dataRef;

            double xMean = 0;
            double yMean = 0;

            // first pass: read in data, compute xbar and ybar
            for (int i=0; i<nPixels; i++) {
                xMean  += x[i] / nPixels;
                yMean  += y[i] / nPixels;
            }
            double xbar = xMean / nPixels;
            double ybar = yMean / nPixels;

            // second pass: compute summary statistics
            double xxbar = 0.0, yybar = 0.0, xybar = 0.0;
            for (int n = 0; n < nPixels; n++) {
                xxbar += (x[n] - xMean) * (x[n] - xMean);
                yybar += (y[n] - yMean) * (y[n] - yMean);
                xybar += (x[n] - xMean) * (y[n] - yMean);
            }

            g = (float) (xybar / xxbar);
            o = (float) (ybar - g * xbar);

            // print results
            //System.out.println("y = " + g + " * x + " + o);

            // normalise data and calculate error map
            for (int n = 0; n < nPixels; n++) {
                dataSRRF[n] = dataSRRF[n] * g + o;
                dataSRRF_CV[n] = dataSRRF_CV[n] * g + o;
                errorMap[n] = abs(dataRef[n] - dataSRRF_CV[n]);
                errorMax = max(errorMap[n], errorMax);
            }
        }
    }
}
