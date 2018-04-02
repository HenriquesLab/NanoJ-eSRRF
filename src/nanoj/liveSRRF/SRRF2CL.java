package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJProfiler;
import nanoj.core2.NanoJThreadExecutor;

import java.io.IOException;
import java.nio.FloatBuffer;
import java.util.ArrayList;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static nanoj.core2.NanoJCL.*;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;

public class SRRF2CL {


    public boolean DEBUG = false;
    private NanoJProfiler prof = new NanoJProfiler();

    private final int nGPUReconstructions;
    public final static int BIN_1 = 1;
    public final static int BIN_2 = 2;
    public final static int BIN_4 = 4;

    public ArrayList<String> reconstructionLabel = new ArrayList<String>();

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

    public ImageStack imsSRRF, imsErrorMap;

    private static int getNGPUReconstructions(int nTimeLags, int doBinFlag) {
        int nReconstructions = 3 + nTimeLags;
        if (doBinFlag == BIN_2) nReconstructions += nTimeLags;
        else if (doBinFlag == BIN_4) nReconstructions += 2*nTimeLags;
        return nReconstructions;
    }

    // return predicted memory usage in MB
    public static double predictMemoryUsed(int width, int height, int nFrames, int magnification, int nTimeLags, int doBinFlag) {
        int nGPUReconstructions = getNGPUReconstructions(nTimeLags, doBinFlag);

        double memUsed = 0;
        memUsed += width * height * nFrames; // clBufferPx
        memUsed += nFrames; // clBufferShiftX
        memUsed += nFrames; // clBufferShiftY
        memUsed += width * height * nFrames; // clBufferGx
        memUsed += width * height * nFrames; // clBufferGy
        memUsed += width * height * magnification * magnification * nGPUReconstructions; // clBufferSRRF
        memUsed += width * height * magnification * magnification; // clBufferSRRF_CVH
        memUsed += width * height * magnification * magnification; // clBufferSRRF_CV

        memUsed *= Float.SIZE / 8000000d;
        return memUsed;
    }

    public SRRF2CL(int width, int height, int nFrames, int magnification, float fwhm, int nTimeLags, int doBinFlag) {
        nGPUReconstructions = getNGPUReconstructions(nTimeLags, doBinFlag);

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
            programString = replaceFirst(programString,"$NTIMELAGS$", ""+nTimeLags);
            programString = replaceFirst(programString,"$DOBIN2$", ""+(doBinFlag == BIN_2 || doBinFlag == BIN_4));
            programString = replaceFirst(programString,"$DOBIN4$", ""+(doBinFlag == BIN_4));
            programString = replaceFirst(programString,"$FWHM$", ""+fwhm);
            programString = replaceFirst(programString,"$SIGMA$", ""+sigma);
            programString = replaceFirst(programString,"$RADIUS$", ""+((int) (sigma * 2) + 1));
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
        clBufferSRRF_CVH = context.createFloatBuffer(widthM * heightM, READ_WRITE);
        clBufferSRRF_CV  = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);

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

        reconstructionLabel.clear();
        reconstructionLabel.add("Raw Average");
        reconstructionLabel.add("SRRF2 Maximum");
        reconstructionLabel.add("SRRF2 Average");

        for (int n=0; n<nTimeLags; n++) reconstructionLabel.add("SRRF2 Accum. TL"+n);
        if (doBinFlag == BIN_2 || doBinFlag == BIN_4)
            for (int n=0; n<nTimeLags; n++) reconstructionLabel.add("SRRF2 Accum. TL"+n+" Bin2");
        if (doBinFlag == BIN_4)
            for (int n=0; n<nTimeLags; n++) reconstructionLabel.add("SRRF2 Accum. TL"+n+" Bin4");

        reconstructionLabel.add("SRRF2 Fusion");
    }

    public synchronized ImageStack calculateSRRF(ImageStack ims, float[] shiftX, float[] shiftY) {
        assert (ims.getWidth() == width && ims.getHeight() == height);
        assert (ims.getSize() <= nFrames);

        int nFrames = ims.getSize();
        int argn;

        float[][] dataSRRF          = new float[nGPUReconstructions][whM];
        float[][] dataSRRFConvolved = new float[nGPUReconstructions][whM];
        float[][] dataErrorMap      = new float[nGPUReconstructions][whM];

        imsSRRF     = new ImageStack(widthM, heightM);
        imsErrorMap = new ImageStack(widthM, heightM);

        // prepare and upload CL buffers
        showStatus("Uploading data to GPU...");
        int id = prof.startTimer();
        fillBuffer(clBufferPx, ims);
        fillBuffer(clBufferShiftX, shiftX);
        fillBuffer(clBufferShiftY, shiftY);
        queue.putWriteBuffer( clBufferPx, false );
        queue.putWriteBuffer( clBufferShiftX, false );
        queue.putWriteBuffer( clBufferShiftY, false );
        prof.recordTime("uploading data to GPU", prof.endTimer(id));

        // make kernelCalculateGradient assignment
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
        // grab frames
        queue.putReadBuffer(clBufferSRRF, true);
        FloatBuffer bufferSRRF = clBufferSRRF.getBuffer();
        for (int r=0; r<nGPUReconstructions; r++) {
            int fOffset = r * whM;
            for (int n=0; n<whM; n++) {
                dataSRRF[r][n] = bufferSRRF.get(fOffset+n);
                if (Float.isNaN(dataSRRF[r][n])) dataSRRF[r][n] = 0; // make sure we dont get any weirdness
            }
        }
        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));
        
        showStatus("Convolving SRRF...");
        float sigmaM = magnification * fwhm / 2.354f;
        id = prof.startTimer();
        for (int r=1; r<nGPUReconstructions; r++) { // note we are ignoring the 0th frame, which is just the interpolated diffraction average
            argn = 0;
            kernelConvolveH.setArg(argn++, clBufferSRRF); // make sure type is the same !!
            kernelConvolveH.setArg(argn++, clBufferSRRF_CVH); // make sure type is the same !!
            kernelConvolveH.setArg(argn++, r); // make sure type is the same !!
            kernelConvolveH.setArg(argn++, sigmaM); // make sure type is the same !!
            queue.put2DRangeKernel(kernelConvolveH, 0, 0, widthM, heightM, 0, 0);
            argn = 0;
            kernelConvolveV.setArg(argn++, clBufferSRRF_CVH); // make sure type is the same !!
            kernelConvolveV.setArg(argn++, clBufferSRRF_CV); // make sure type is the same !!
            kernelConvolveV.setArg(argn++, sigmaM); // make sure type is the same !!
            queue.put2DRangeKernel(kernelConvolveV, 0, 0, widthM, heightM, 0, 0);
            // grab convolved frame
            queue.putReadBuffer(clBufferSRRF_CV, true);
            grabBuffer(clBufferSRRF_CV, dataSRRFConvolved[r], true);
        }
        prof.recordTime("kernelConvolveH & kernelConvolveV", prof.endTimer(id));

        id = prof.startTimer();
        queue.finish(); // make sure everything is done...
        prof.recordTime("waiting for queue to finish", prof.endTimer(id));

        // grab interpolated data frame
        showStatus("Preparing data...");
        float[] RAW_AVE = dataSRRF[0];

        // grab data from convolved frames
        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        for (int r=1; r<nGPUReconstructions; r++) {
            // this thread will read the buffer for SRRF and SRRF_CV plus rescale their intensity to match RAW_AVE
            ThreadedNormaliseAndCalculateErrorMap t = new ThreadedNormaliseAndCalculateErrorMap(RAW_AVE, dataSRRF[r], dataSRRFConvolved[r], dataErrorMap[r]);
            NTE.execute(t);
        }
        NTE.finish();

        // calculate the error map and send off each reconstruction to imsSRRF
        for (int r=0; r<nGPUReconstructions; r++) {
            imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF[r]));
            imsErrorMap.addSlice(new FloatProcessor(widthM, heightM, dataErrorMap[r]));
        }

        showStatus("Calculating SRRF Fusion...");
        float[] fusion = new float[whM];
        double[] weightSum = new double[whM];

        for (int n = 0; n < whM; n++) {
            float errorMax = -Float.MAX_VALUE;
            for (int r=1; r<nGPUReconstructions; r++) errorMax = max(dataErrorMap[r][n], errorMax);

            for (int r=1; r<nGPUReconstructions; r++) {
                float v = dataSRRF[r][n];
                double w = errorMax / max(dataErrorMap[r][n], 1);
                //double w = 1 / max(errorMap[r][n], 1);
                fusion[n] += v * w;
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

    class ThreadedNormaliseAndCalculateErrorMap extends Thread {
        private final float[] dataRef;
        private final float[] dataSRRF;
        private final float[] dataSRRFConvolved;
        private final float[] dataErrorMap;

        public float g, o;

        public ThreadedNormaliseAndCalculateErrorMap(float[] dataRef, float dataSRRF[], float[] dataSRRFConvolved, float[] dataErrorMap) {
            this.dataRef = dataRef;
            this.dataSRRF = dataSRRF;
            this.dataSRRFConvolved = dataSRRFConvolved;
            this.dataErrorMap = dataErrorMap;
        }

        @Override
        public void run() {
            int nPixels = dataRef.length;

            // now calculate linear regression
            float[] x = dataSRRFConvolved;
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
                xybar += (x[n] - xMean) * (y[n] - yMean);
            }

            g = (float) (xybar / xxbar);
            o = (float) (ybar - g * xbar);

            // normalise data and calculate error map
            for (int n = 0; n < nPixels; n++) {
                dataSRRF[n] = dataSRRF[n] * g + o;
                dataSRRFConvolved[n] = dataSRRFConvolved[n] * g + o;
                dataErrorMap[n] = abs(dataRef[n] - dataSRRFConvolved[n]);
            }
        }
    }
}
