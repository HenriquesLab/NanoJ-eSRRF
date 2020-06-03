package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.plugin.ZAxisProfiler;
import ij.process.ImageProcessor;
import nanoj.liveSRRF.SSIMCalculator;
import nanoj.core2.NanoJPrefs;


import java.util.Arrays;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;


public class SSIMtrace_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private ZAxisProfiler zProfiler = new ZAxisProfiler();


    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();


        IJ.log("\\Clear");  // Clear the log window

        ImageStack ims = imp.getStack();

        int nFrames = imp.getImageStack().getSize();
        int width = imp.getImageStack().getWidth();
        int height = imp.getImageStack().getHeight();
        String imageTitle = imp.getTitle();

        String[] regularizationMethods = new String[2];
        int n = 0;
        regularizationMethods[n++] = "image bit depth";
        regularizationMethods[n++] = "image dynamic range";

        // ----- GUI -----
        GenericDialog gd = new GenericDialog("liveSRRF - SSIM estimator");
        gd.addChoice("Regularization method", regularizationMethods, prefs.get("regMethod", regularizationMethods[0]));
        gd.addNumericField("Radius (in pixels)", prefs.get("radius", 2), 2);
        gd.addNumericField("Frames to average for ref.", prefs.get("nFramesToAverage", 5), 0);
        gd.addCheckbox("Temporal blocking", prefs.get("tempBlocking", true));
        gd.addNumericField("Block size (frames)", prefs.get("blockSize", 100), 0);
        gd.addCheckbox("Show trace", prefs.get("showTrace", true));
        gd.addNumericField("Smoothing factor", prefs.get("smoothingFactor", 0.9f),2);
        gd.addNumericField("Sigma cut-off", prefs.get("sigmaCutOff", 10),1);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String regMethod = gd.getNextChoice();
        double radius = gd.getNextNumber();
        int nFramesToAverage = (int) gd.getNextNumber();
        boolean tempBlocking = gd.getNextBoolean();
        int blockSize = (int) gd.getNextNumber();
        boolean showTrace = gd.getNextBoolean();
        float smoothingFactor = (float) gd.getNextNumber();
        float sigmaCutOff = (float) gd.getNextNumber();

        prefs.set("regMethod", regMethod);
        prefs.set("radius", (float) radius);
        prefs.set("nFramesToAverage", nFramesToAverage);
        prefs.set("tempBlocking", tempBlocking);
        prefs.set("blockSize", blockSize);
        prefs.set("showTrace", showTrace);
        prefs.set("smoothingFactor", smoothingFactor);
        prefs.set("sigmaCutOff", sigmaCutOff);

        int nBlocks = 1;
        if (tempBlocking) {
            nBlocks = (int) (float) nFrames / blockSize;
        }
        else{
            blockSize = nFrames;
        }


        // Get the regularization factors from the stack
        float[] regFactors = getRegularizationFactors(ims, regMethod);
        int startingPoint, currentBlockSize;
        float[] ssimTrace;
        double[] frames;
        double[][] dataSmoothed;
        int cutOffPosition;

        for (int b = 0; b < nBlocks; b++) {

            IJ.log("Block #"+(b+1));
            // Get the reference image by averaging the first "nFramesToAverage" frames
            ImageProcessor ipRef = ims.duplicate().getProcessor(1 + b*blockSize).convertToFloatProcessor();
            for (int i = 1; i < nFramesToAverage; i++) {
                ipRef = add(ipRef.convertToFloatProcessor(), ims.duplicate().getProcessor(i + 1 + b*blockSize).convertToFloatProcessor());
            }
            ipRef.multiply(1.0f / (double) nFramesToAverage);


            SSIMCalculator ssimCalculator = new SSIMCalculator(ipRef, radius, regFactors);
            ImageStack imsSSIM = new ImageStack(width, height, blockSize);
            startingPoint = 1;
            currentBlockSize = blockSize;
            if (b == 0){
                startingPoint = nFramesToAverage+1;
                imsSSIM = new ImageStack(width, height, blockSize-nFramesToAverage);
                currentBlockSize = blockSize - nFramesToAverage;
            }
            for (int i = 0; i < currentBlockSize; i++) {
//                IJ.showProgress(i, currentBlockSize);
                ImageProcessor ip = ims.getProcessor(startingPoint + i + b*blockSize);
                imsSSIM.setProcessor(ssimCalculator.Calculate(ip), i+1);
            }

            ImagePlus impSSIM = new ImagePlus();
            impSSIM.setStack(imageTitle + " - SSIM stack (block #"+(b+1)+")", imsSSIM);
            impSSIM.show();


            Plot ssimPlot = zProfiler.getPlot(impSSIM);
            ssimTrace = ssimPlot.getYValues();

            frames = new double[ssimTrace.length];
            for (int f = 0; f < ssimTrace.length; f++) {
                frames[f] = f+1;
            }

            dataSmoothed = exponentialCurveSmoothing(ssimTrace, smoothingFactor, true);
            cutOffPosition = getSigmacutOffPosition(ssimTrace, smoothingFactor, sigmaCutOff);

            if (showTrace) {
                ssimPlot.setAxisXLog(true);
                ssimPlot.setXYLabels("Number of frames", "SSIM score");
//                ssimPlot.setLimits(NaN, NaN, -0.1, 1.0);
                ssimPlot.show();
                ssimPlot.setColor("red");
                ssimPlot.setLineWidth(2);
                ssimPlot.add("line", frames, dataSmoothed[0]);
                ssimPlot.setColor("blue");
                ssimPlot.add("line", new double[]{cutOffPosition, cutOffPosition}, new double[]{-1, 1});

                Plot deltaPlot = new Plot(imageTitle + " - SSIM stack (block #"+(b+1)+") Delta Plot", "Number of frames", "Difference of SSIM score");
                deltaPlot.add("line", frames, dataSmoothed[1]);
                deltaPlot.setAxisXLog(true);
                deltaPlot.show();
            }


        }
    }

    // ------------------ Functions ------------------
    private static float[] getRegularizationFactors(ImageStack ims, String regMethod){

        float[] C = new float[2];
        // Regularization factors
        if (regMethod.equals("image bit depth")) {
            int bitDepth = ims.getBitDepth();
            C[0] = (float) Math.pow(0.01f * (float) (Math.pow(2.0f, bitDepth) - 1), 2);
            IJ.log("Bit depth: "+bitDepth);
        }
        else{
            float[] percentiles = new float[2];
            percentiles[0] = getPercentileFromImageStack(ims, 0.01f);
            percentiles[1] = getPercentileFromImageStack(ims, 0.99f);
            IJ.log("1% percentile: "+percentiles[0]);
            IJ.log("99% percentile: "+percentiles[1]);
            C[0] = (float) Math.pow(0.01f * (percentiles[1] - percentiles[0]), 2);
        }
        C[1] = 9 * C[0];

        IJ.log("Regularization factors: "+regMethod);
        IJ.log("C1: "+C[0]);
        IJ.log("C2: "+C[1]);

        return C;

    }

    public static float getPercentileFromImageStack(ImageStack ims, float percentile){

        float[] pixelValues = ImageStackToFloatArray(ims);
        Arrays.sort(pixelValues);
        int nPixels = ims.getWidth()*ims.getHeight()*ims.getSize();
        int id = (int) ((float) nPixels*percentile);

        return pixelValues[id];
    }

    public static double[][] exponentialCurveSmoothing(float[] data, float smoothingfactor, boolean flip){

        double[][] data_smoothed = new double[2][data.length];

        if (flip) {
            data_smoothed[0][data.length-1] = data[data.length-1];
            data_smoothed[1][data.length-1] = data[data.length-1] - data_smoothed[0][data.length-1];
            for (int i = 1; i < data.length; i++) {
                data_smoothed[0][data.length-1 - i] = smoothingfactor * data[data.length-1 - i] + (1.0f - smoothingfactor) * data_smoothed[0][data.length-1 - (i-1)];
                data_smoothed[1][data.length-1 - i] = data[data.length-1 - i] - data_smoothed[0][data.length-1 - i];
            }
        }
        else{
            data_smoothed[0][0] = data[0];
//            data_smoothed[0][0] = smoothingfactor*data[0] + (1.0f-smoothingfactor)*1.0f; // assuming previous data was 1
            data_smoothed[1][0] = data[0] - data_smoothed[0][0];
            for (int i = 1; i < data.length; i++) {
                data_smoothed[0][i] = smoothingfactor * data[i] + (1.0f - smoothingfactor) * data_smoothed[0][i - 1];
                data_smoothed[1][i] = data[i] - data_smoothed[0][i];
            }
        }

        return data_smoothed;
    }

    public static int getSigmacutOffPosition(float[] data, float smoothingfactor, float sigmaCutOff){

        double[][] data_smoothed = exponentialCurveSmoothing(data, smoothingfactor, true);

        float m = 0, s = 0;
        for (int i = 0; i < data_smoothed[0].length; i++) {
            m += data_smoothed[1][i];
            s += data_smoothed[1][i]*data_smoothed[1][i];
        }
        float std = (float) Math.sqrt(s/data_smoothed[0].length - (m/data_smoothed[0].length)*(m/data_smoothed[0].length));

        int cutOffPosition = 0;
        while ((data_smoothed[0][cutOffPosition] > data_smoothed[0][0]-sigmaCutOff*std) && (cutOffPosition<data_smoothed[0].length-1)) {
            cutOffPosition ++;
        }
        cutOffPosition ++; // adjust for the fact that frame numbers start at 1
        IJ.log("Sigma: "+std);
        IJ.log("Cut-off position: "+cutOffPosition);
        return cutOffPosition;
    }

}