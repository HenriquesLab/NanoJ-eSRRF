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
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.*;

public class SRRF2_CL {

    private NanoJProfiler prof = new NanoJProfiler();

    public ArrayList<String> reconstructionLabel = new ArrayList<String>();

    static CLContext context;
    static CLProgram programLiveSRRF, programConvolve2DIntegratedGaussian;
    static CLKernel kernelCalculateGradient, kernelCalculateSRRF, kernelConvolveH, kernelConvolveV;
    static CLCommandQueue queue;

    private final int width, height, widthM, heightM, nFrameOnGPU, magnification, whM;
    private final float fwhm;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferShiftX, clBufferShiftY,
            clBufferSRRF;

    public ImageStack imsSRRF, imsErrorMap;


    public SRRF2_CL(int width, int height, int nFramesOnGPU, int magnification, float fwhm, int sensitivity) {

        this.nFrameOnGPU = nFramesOnGPU;
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
            float sigma = fwhm / 2.354f;
            String programString = getResourceAsString(SRRF2_CL.class, "liveSRRF.cl");
            programString = replaceFirst(programString, "$MAX_FRAMES$", "" + nFramesOnGPU);
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

            programConvolve2DIntegratedGaussian = context.createProgram(SRRF2_CL.class.getResourceAsStream("/Convolve2DIntegratedGaussian.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }
        kernelCalculateGradient = programLiveSRRF.createCLKernel("calculateGradient_2point");
        kernelCalculateSRRF = programLiveSRRF.createCLKernel("calculateSRRF");
        kernelConvolveH = programConvolve2DIntegratedGaussian.createCLKernel("convolveHorizontal");
        kernelConvolveV = programConvolve2DIntegratedGaussian.createCLKernel("convolveVertical");

        clBufferPx = context.createFloatBuffer(width * height * nFramesOnGPU, READ_ONLY);
        clBufferShiftX = context.createFloatBuffer(nFramesOnGPU, READ_ONLY);
        clBufferShiftY = context.createFloatBuffer(nFramesOnGPU, READ_ONLY);
        clBufferGx = context.createFloatBuffer(4 * width * height * nFramesOnGPU, READ_WRITE);
        clBufferGy = context.createFloatBuffer(4 * width * height * nFramesOnGPU, READ_WRITE);
        clBufferSRRF = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);

        System.out.println("used device memory: " + (
                clBufferPx.getCLSize() +
                        clBufferShiftX.getCLSize() +
                        clBufferShiftY.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferSRRF.getCLSize())
                / 1000000d + "MB");

        reconstructionLabel.clear();
        reconstructionLabel.add("Raw Average");
        reconstructionLabel.add("liveSRRF Average");
        reconstructionLabel.add("liveSRRF Fusion");
    }

    public synchronized ImageStack calculateSRRF(ImageStack ims, float[] shiftX, float[] shiftY) {
        assert (ims.getWidth() == width && ims.getHeight() == height);
        assert (ims.getSize() <= nFrameOnGPU);

        int nFrames = ims.getSize();
        int argn;

        float[] dataSRRF = new float[whM];
        float[] dataSRRFConvolved = new float[whM];
        float[] dataErrorMap = new float[whM];

        imsSRRF = new ImageStack(widthM, heightM);
        imsErrorMap = new ImageStack(widthM, heightM);

        // prepare and upload CL buffers
        showStatus("Uploading data to GPU...");
        int id = prof.startTimer();
        fillBuffer(clBufferPx, ims);
        fillBuffer(clBufferShiftX, shiftX);
        fillBuffer(clBufferShiftY, shiftY);
        queue.putWriteBuffer(clBufferPx, false);
        queue.putWriteBuffer(clBufferShiftX, false);
        queue.putWriteBuffer(clBufferShiftY, false);
        prof.recordTime("uploading data to GPU", prof.endTimer(id));

        // make kernelCalculateGradient assignment
        argn = 0;
        showStatus("Calculating gradient...");
        id = prof.startTimer();
        kernelCalculateGradient.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelCalculateGradient.setArg(argn++, nFrames); // make sure type is the same !!
        queue.put2DRangeKernel(kernelCalculateGradient, 0, 0, width, height, 0, 0);
        prof.recordTime("kernelCalculateGradient", prof.endTimer(id));

        // make kernelCalculateSRRF assignment
        argn = 0;
        showStatus("Calculating SRRF...");
        id = prof.startTimer();
        kernelCalculateSRRF.setArg(argn++, clBufferPx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGx); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferGy); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftX); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferShiftY); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, clBufferSRRF); // make sure type is the same !!
        kernelCalculateSRRF.setArg(argn++, nFrames); // make sure type is the same !!
        // break calculation into 128x128 blocks to make them more stable
        int nXBlocks = widthM / 128 + ((widthM % 128 == 0) ? 0 : 1);
        int nYBlocks = heightM / 128 + ((heightM % 128 == 0) ? 0 : 1);
        for (int nYB = 0; nYB < nYBlocks; nYB++) {
            int yWorkSize = min(128, heightM - nYB * 128);
            for (int nXB = 0; nXB < nXBlocks; nXB++) {
                showStatus("Calculating SRRF... blockX=" + nXB + "/" + nXBlocks + "  blockY=" + nYB + "/" + nYBlocks);
                int xWorkSize = min(128, widthM - nXB * 128);
                queue.put2DRangeKernel(kernelCalculateSRRF, nXB * 128, nYB * 128, xWorkSize, yWorkSize, 0, 0);
            }
        }
        // grab frames
        queue.putReadBuffer(clBufferSRRF, true);
        FloatBuffer bufferSRRF = clBufferSRRF.getBuffer();

        for (int n = 0; n < whM; n++) {
            dataSRRF[n] = bufferSRRF.get(whM + n);
            if (Float.isNaN(dataSRRF[n])) dataSRRF[n] = 0; // make sure we don't get any weirdness
        }

        prof.recordTime("kernelCalculateSRRF", prof.endTimer(id));


        id = prof.startTimer();
        queue.finish(); // make sure everything is done...
        prof.recordTime("waiting for queue to finish", prof.endTimer(id));

        // grab interpolated data frame
        showStatus("Preparing data...");
        float[] RAW_AVE = dataSRRF;

        // grab data from convolved frames
        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        // this thread will read the buffer for SRRF and SRRF_CV plus rescale their intensity to match RAW_AVE
        ThreadedNormaliseAndCalculateErrorMap t = new ThreadedNormaliseAndCalculateErrorMap(RAW_AVE, dataSRRF, dataSRRFConvolved, dataErrorMap, widthM, heightM);
        NTE.execute(t);

        NTE.finish();

        // calculate the error map and send off each reconstruction to imsSRRF
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, dataSRRF));
        imsErrorMap.addSlice(new FloatProcessor(widthM, heightM, dataErrorMap));


        showStatus("Calculating SRRF Fusion...");
        float[] fusion = new float[whM];
        double[] weightSum = new double[whM];

        for (int n = 0; n < whM; n++) {
            float errorMax = -Float.MAX_VALUE;
            errorMax = max(dataErrorMap[n], errorMax);

            float v = dataSRRF[n];
            double w = errorMax / max(dataErrorMap[n], 1);
            //double w = 1 / max(errorMap[r][n], 1);
            fusion[n] += v * w;
            weightSum[n] += w;
        }

        for (int n = 0; n < whM; n++) fusion[n] /= weightSum[n];
        imsSRRF.addSlice(new FloatProcessor(widthM, heightM, fusion));

        return imsSRRF;
    }

    public void showStatus(String text) {
        IJ.showStatus("liveSRRF: " + text);
    }

    public void release() {
        context.release();
    }

    class ThreadedNormaliseAndCalculateErrorMap extends Thread {
        private final float[] dataRef;
        private final float[] dataSRRF;
        private final float[] dataSRRFConvolved;
        private final float[] dataErrorMap;
        private final int width;
        private final int height;

        public float g, o;

        public ThreadedNormaliseAndCalculateErrorMap(float[] dataRef, float dataSRRF[], float[] dataSRRFConvolved, float[] dataErrorMap, int width, int height) {
            this.dataRef = dataRef;
            this.dataSRRF = dataSRRF;
            this.dataSRRFConvolved = dataSRRFConvolved;
            this.dataErrorMap = dataErrorMap;
            this.width = width;
            this.height = height;
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
            for (int i = 0; i < nPixels; i++) {
                xMean += x[i] / nPixels;
                yMean += y[i] / nPixels;
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

            float[] dataErrorMapSmooth = dataErrorMap.clone();

            // smooth error map, to minimize corruption by noise
            for (int j = 1; j < height - 1; j++) {
                int l0 = (j - 1) * width;
                int l1 = j * width;
                int l2 = (j + 1) * width;

                for (int i = 1; i < width - 1; i++) {
                    float vMean = 0;
                    vMean += dataErrorMap[l0 + i - 1];
                    vMean += dataErrorMap[l0 + i];
                    vMean += dataErrorMap[l0 + i + 1];
                    vMean += dataErrorMap[l1 + i - 1];
                    vMean += dataErrorMap[l1 + i];
                    vMean += dataErrorMap[l1 + i + 1];
                    vMean += dataErrorMap[l2 + i - 1];
                    vMean += dataErrorMap[l2 + i];
                    vMean += dataErrorMap[l2 + i + 1];
                    vMean /= 9;
                    dataErrorMapSmooth[l1 + i] = vMean;
                }
            }
            System.arraycopy(dataErrorMapSmooth, 0, dataErrorMap, 0, dataErrorMap.length);
        }
    }


}
