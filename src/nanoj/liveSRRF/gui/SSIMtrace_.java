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

import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.liveSRRF.SSIMCalculator.getRegularizationFactors;
import static org.apache.commons.math.stat.inference.TestUtils.tTest;


public class SSIMtrace_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
//    private ZAxisProfiler zProfiler = new ZAxisProfiler();


    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();


        IJ.log("\\Clear");  // Clear the log window

        ImageStack ims = imp.getStack();

        int nFrames = imp.getImageStack().getSize();
        String imageTitle = imp.getTitle();

        String[] regularizationMethods = new String[2];
        int n = 0;
        regularizationMethods[n++] = "image bit depth";
        regularizationMethods[n++] = "image dynamic range";

        // ----- GUI -----
        GenericDialog gd = new GenericDialog("eSRRF - SSIM estimator");
        gd.addChoice("Regularization method", regularizationMethods, prefs.get("regMethod", regularizationMethods[0]));
//        gd.addNumericField("Radius (in pixels)", prefs.get("radius", 2), 2);
        gd.addNumericField("Block size (frames, 0 = auto)", prefs.get("blockSize", 200), 0);
//        gd.addCheckbox("Show trace", prefs.get("showTrace", true));
        gd.addNumericField("Smoothing factor", prefs.get("smoothingFactor", 0.15f),2);
        gd.addNumericField("Sigma cut-off", prefs.get("sigmaCutOff", 2),1);
        gd.addCheckbox("Calculate Cut-off over time", prefs.get("calculateCutoffVStime", true));
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String regMethod = gd.getNextChoice();
//        double radius = gd.getNextNumber();
        int blockSize = (int) gd.getNextNumber();
//        boolean showTrace = gd.getNextBoolean();
        float smoothingFactor = (float) gd.getNextNumber();
        float sigmaCutOff = (float) gd.getNextNumber();
        boolean calculateCutoffVStime = gd.getNextBoolean();

        prefs.set("regMethod", regMethod);
//        prefs.set("radius", (float) radius);
        prefs.set("blockSize", blockSize);
//        prefs.set("showTrace", showTrace);
        prefs.set("smoothingFactor", smoothingFactor);
        prefs.set("sigmaCutOff", sigmaCutOff);
        prefs.set("calculateCutoffVStime", calculateCutoffVStime);



        // Get the regularization factors from the stack
        float[] regFactors = getRegularizationFactors(ims, regMethod);
        SSIMCalculator ssimCalculator;

        int nBlocks = nFrames - blockSize - 1;

        if  (blockSize == 0) {
            ssimCalculator = new SSIMCalculator(ims.duplicate().getProcessor(1).convertToFloatProcessor(), 0, regFactors);

            double[] initialSSIMtrace = new double[nFrames - 1];
            for (int i = 0; i < nFrames - 1; i++) {
                initialSSIMtrace[i] = ssimCalculator.CalculateMetric(ims.getProcessor(i + 2).convertToFloatProcessor());
            }

            int initialCutOffPosition = getSigmacutOffPosition(initialSSIMtrace, smoothingFactor, sigmaCutOff, true);
            IJ.log("Initial Cut-off position: " + initialCutOffPosition);
            blockSize = Math.min(2*initialCutOffPosition, nFrames-1);
            nBlocks = nFrames - blockSize - 1;
        }

        if (nBlocks < 0) {
            nBlocks = 1;
            blockSize = nFrames-1;
        }

        IJ.log("Block size: "+blockSize);
        IJ.log("Number of blocks: "+nBlocks);
        IJ.log("Sigma cut-off: "+sigmaCutOff);

//        ImageProcessor ipRef;
        double[][] ssimTraces = new double[nBlocks][blockSize];

        for (int b = 0; b < nBlocks; b++) {
            IJ.showProgress(b, nBlocks);

//            ipRef = ims.duplicate().getProcessor(1 + b).convertToFloatProcessor(); // rolling analysis
            ssimCalculator = new SSIMCalculator(ims.getProcessor(1 + b).convertToFloatProcessor(), Float.MAX_VALUE, regFactors);

//            ImageStack imsBlock = new ImageStack(ims.getWidth(), ims.getHeight());
//            for (int i = 0; i < blockSize; i++) {
//                imsBlock.addSlice(ims.getProcessor(b+i+2).convertToFloatProcessor());
//            }
//            ssimTraces[b] = ssimCalculator.CalculateMetric(imsBlock);

            for (int i = 0; i < blockSize; i++) {
                ssimTraces[b][i] = ssimCalculator.CalculateMetric(ims.getProcessor(b+i+2).convertToFloatProcessor());
            }
        }

        displaySSIMtraces(ssimTraces, imageTitle, sigmaCutOff);
        if (calculateCutoffVStime) getTimeConstantOverTime(ssimTraces, smoothingFactor, sigmaCutOff);

    }

    // ------------------ Functions ------------------

    public static double[][] exponentialCurveSmoothing(double[] data, float smoothingfactor, boolean flip){

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

    public static int getSigmacutOffPosition(double[] data, float smoothingfactor, float sigmaCutOff, boolean showTrace){

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

        if (showTrace) {
            IJ.log("Sigma: "+std);
            double[] frames = new double[data.length];
            for (int f = 0; f < data.length; f++) {
                frames[f] = f+1;
            }
            Plot tracePlot = new Plot( "SSIM plot smoothed curve", "Number of frames", "Difference of SSIM score");
            tracePlot.add("line", frames, data);
            tracePlot.setAxisXLog(true);

//                ssimPlot.setLimits(NaN, NaN, -0.1, 1.0);
            tracePlot.show();
            tracePlot.setColor("red");
            tracePlot.setLineWidth(2);
            tracePlot.add("line", frames, data_smoothed[0]);
            tracePlot.setColor("blue");
            tracePlot.add("line", new double[]{cutOffPosition, cutOffPosition}, new double[]{-1, 1});
            tracePlot.show();
        }

        return cutOffPosition;
    }


    public static int getSigmacutOffPosition(double[] meanArray, double[] stdArray, float sigmaCutOff){

        double meanSTD = 0;
        for (int i = 0; i < meanArray.length; i++) {
            meanSTD += stdArray[i];
        }
        meanSTD /= meanArray.length;

        int cutOffPosition = 0;
        while ((meanArray[cutOffPosition] > meanArray[0]-sigmaCutOff*meanSTD) && (cutOffPosition<meanArray.length-1)) {
            cutOffPosition ++;
        }
        cutOffPosition ++; // adjust for the fact that frame numbers start at 1
        IJ.log("Sigma: "+meanSTD);
        IJ.log("Cut-off position: "+cutOffPosition);
        return cutOffPosition;
    }


    public static void displaySSIMtraces(double[][] ssimTraces, String title, float sigmaCutOff){

        Plot ssimTracePlot = new Plot(title + " - SSIM metric", "Number of frames", "SSIM metric");
        int traceLength = ssimTraces[0].length;
        int nTraces = ssimTraces.length;

        double[] frames = new double[traceLength];
        for (int f = 0; f < traceLength; f++) {
            frames[f] = f+1;
        }

        for (int i = 0; i < nTraces; i++) {
            ssimTracePlot.add("line", frames, ssimTraces[i]);
        }

        ssimTracePlot.setAxisXLog(true);
        ssimTracePlot.show();

        double[][] meanSTDptrace = getMeanSTDPvaluefromArray(ssimTraces);
        double[][] stdOffsetTraces = new double[2][traceLength];
        for (int i = 0; i < traceLength; i++) {
            stdOffsetTraces[0][i] = meanSTDptrace[0][i] + meanSTDptrace[1][i];
            stdOffsetTraces[1][i] = meanSTDptrace[0][i] - meanSTDptrace[1][i];
        }

        int cutoffPosition = getSigmacutOffPosition(meanSTDptrace[0], meanSTDptrace[1], sigmaCutOff);

        ssimTracePlot.setColor("red");
        ssimTracePlot.setLineWidth(2);
        ssimTracePlot.add("line", frames, meanSTDptrace[0]);
        ssimTracePlot.setLineWidth(1);
        ssimTracePlot.add("line", frames, stdOffsetTraces[0]);
        ssimTracePlot.add("line", frames, stdOffsetTraces[1]);
        ssimTracePlot.setColor("blue");
        ssimTracePlot.add("line", new double[]{cutoffPosition, cutoffPosition}, new double[]{-1, 1});
        ssimTracePlot.show();
    }


    public static double[][] getMeanSTDPvaluefromArray(double[][] traces){
        int traceLength = traces[0].length;
        int nTraces = traces.length;

        double[][] meanSTDPvalTraces = new double[3][traceLength];

//        double[] sampleT0 = new double[nTraces];
//        for (int j = 0; j < nTraces; j++) {
//            sampleT0[j] = traces[j][0];
//        }

//        double[] sampleT = new double[nTraces];
        for (int i = 0; i < traceLength; i++) {
            for (int j = 0; j < nTraces; j++) {
                meanSTDPvalTraces[0][i] += traces[j][i]/nTraces;
                meanSTDPvalTraces[1][i] += traces[j][i]*traces[j][i]/nTraces;
//                sampleT[j] = traces[j][i];
            }
//            try{
//                meanSTDPvalTraces[2][i] = Math.log10(tTest(sampleT0, sampleT));}
//            catch(Exception e){
////                IJ.log("T-test unsuccessful...");
//            }
//            if (meanSTDPvalTraces[2][i] == Double.MIN_VALUE){
//                IJ.log("Prout!");
//                meanSTDPvalTraces[2][i] = Double.NaN;
//            }
        }

        for (int i = 0; i < traceLength; i++) {
            meanSTDPvalTraces[1][i] = Math.sqrt(meanSTDPvalTraces[1][i] - meanSTDPvalTraces[0][i]*meanSTDPvalTraces[0][i]);
        }

        return meanSTDPvalTraces;
    }


    public static double[] getTimeConstantOverTime(double[][] data, float smoothingfactor, float sigmaCutOff){

        int nTraces = data.length;
//        int traceLength = data[0].length;
        double[] cutOffList = new double[nTraces];

        for (int i = 0; i < nTraces; i++) {
            cutOffList[i] = getSigmacutOffPosition(data[i], smoothingfactor, sigmaCutOff, false);
        }

        double[] frames = new double[nTraces];
        for (int f = 0; f < nTraces; f++) {
            frames[f] = f+1;
        }

        Plot cutoffPlot = new Plot("Cut-off vs. time", "Frames ", "Time constant (frames)");
        cutoffPlot.add("line", frames, cutOffList);
        cutoffPlot.show();

        return cutOffList;
    }




}