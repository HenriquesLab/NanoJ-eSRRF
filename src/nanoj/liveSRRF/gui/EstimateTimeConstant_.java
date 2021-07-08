package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import nanoj.liveSRRF.OpticalFlowMagnitudeEstimator;
import nanoj.liveSRRF.SSIMCalculator;

import java.util.Arrays;

import static java.lang.Float.isNaN;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;


public class EstimateTimeConstant_ implements PlugIn {

    private final boolean DEBUG = false;

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

        String[] estimationMethods = new String[2];
        int n = 0;
        estimationMethods[n++] = "SSIM";
        estimationMethods[n++] = "Optical flow magnitude";

        String[] regularizationMethods = new String[2];
        n = 0;
        regularizationMethods[n++] = "image bit depth";
        regularizationMethods[n++] = "image dynamic range";

        GenericDialog gd = new GenericDialog("eSRRF - Time constant estimator");
        gd.addChoice("Estimation method", estimationMethods, estimationMethods[1]);

        gd.addMessage("---- SSIM method ----");
        gd.addChoice("Regularization method", regularizationMethods, regularizationMethods[0]);
        gd.addNumericField("Radius (in pixels)", 2, 2);

        gd.addMessage("---- Optical flow method ----");
        gd.addNumericField("Noise standard deviation (ADC)", 65, 1);

        gd.addMessage("---- Main parameters ----");
        gd.addNumericField("Frames to average for ref.", 5, 0);
        gd.addCheckbox("Temporal blocking", false);
        gd.addNumericField("Block size (frames)", 100, 0);
        gd.addNumericField("Percentile for plot: ", 0.90,2);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String estimationMethod = gd.getNextChoice();
        String regMethod = gd.getNextChoice();
        double radius = gd.getNextNumber();
        float noiseStd = (float) gd.getNextNumber();

        int nFramesToAverage = (int) gd.getNextNumber();
        boolean tempBlocking = gd.getNextBoolean();
        int blockSize = (int) gd.getNextNumber();
        float prctilePlot = (float) gd.getNextNumber();

        int nBlocks = 1;
        if (tempBlocking) {
            nBlocks = (int) (float) nFrames / blockSize;
        }
        else{
            blockSize = nFrames;
        }

        // Main loop
        for (int b = 0; b < nBlocks; b++) {
            IJ.log("Block #"+(b+1));
            // Get the reference frame
            ImageProcessor ipRef = ims.duplicate().getProcessor(1 + b*blockSize).convertToFloatProcessor();
            for (int i = 1; i < nFramesToAverage; i++) {
                ipRef = add(ipRef.convertToFloatProcessor(), ims.duplicate().getProcessor(i + 1 + b*blockSize).convertToFloatProcessor());
            }
            ipRef.multiply(1.0f / (double) nFramesToAverage);

            ImageStack imsMotion = new ImageStack(width, height, blockSize);
            if (estimationMethod.equals("SSIM")){
                float[] regFactors = getRegularizationFactors(ims, regMethod);
                SSIMCalculator ssimCalculator = new SSIMCalculator(ipRef, radius, regFactors);
                for (int i = 1; i <= blockSize; i++) {
                    IJ.showProgress(i, blockSize);
                    ImageProcessor ip = ims.getProcessor(i + b*blockSize);
                    imsMotion.setProcessor(ssimCalculator.CalculateSSIMmap(ip), i);
                }
            }

            if (estimationMethod.equals("Optical flow magnitude")){
                OpticalFlowMagnitudeEstimator flowEstimator = new OpticalFlowMagnitudeEstimator(ipRef, noiseStd);
                for (int i = 1; i <= blockSize; i++) {
                    IJ.showProgress(i, blockSize);
                    ImageProcessor ip = ims.getProcessor(i + b*blockSize);
                    imsMotion.setProcessor(flowEstimator.calculateMagnitude(ip), i);
                }
            }

            ImagePlus impMotion = new ImagePlus();
            impMotion.setStack(imageTitle + " - "+estimationMethod+" stack (block #"+(b+1)+")", imsMotion);
            impMotion.show();

            double[] prctileArray = getPercentileFromImageStack(imsMotion, prctilePlot, false);
            Plot thisPlot = new Plot("Motion metrics", "Frame #", "Metric (AU)");
            thisPlot.add("connected circle", prctileArray);
            thisPlot.show();

        }

        IJ.log("-------------------");
        IJ.log("All done.");
    }

    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------

    private static float[] getRegularizationFactors(ImageStack ims, String regMethod){
        float[] C = new float[2];
        // Regularization factors
        if (regMethod.equals("image bit depth")) {
            int bitDepth = ims.getBitDepth();
            C[0] = (float) Math.pow(0.01f * (float) (Math.pow(2.0f, bitDepth) - 1), 2);
            IJ.log("Bit depth: "+bitDepth);
        }
        else{
            double[] percentiles = new double[2];
            percentiles[0] = getPercentileFromImageStack(ims, 0.01f, true)[0];
            percentiles[1] = getPercentileFromImageStack(ims, 0.99f, true)[0];
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

    public static double[] getPercentileFromImageStack(ImageStack ims, float percentile, boolean combine){

        if (combine) {
            float[] pixelValues = ImageStackToFloatArray(ims);
            Arrays.sort(pixelValues);
            int nPixels = ims.getWidth()*ims.getHeight()*ims.getSize();
            int id = (int) ((float) nPixels*percentile);

            return new double[]{pixelValues[id]};
        }
        else {

            double[] prctileArray = new double[ims.getSize()];
            for (int i = 0; i < ims.getSize(); i++) {
                float[] pixelValues = (float[]) ims.getProcessor(i+1).duplicate().convertToFloatProcessor().getPixels();

//                IJ.log("Size: "+pixelValues.length);
                Arrays.sort(pixelValues);
                int nPixels = ims.getWidth()*ims.getHeight();

                int nPixelsIgnoringNaN = nPixels-1;
                while (isNaN(pixelValues[nPixelsIgnoringNaN])) nPixelsIgnoringNaN--;
//                IJ.log("Number of pixels that are not NaN: "+nPixelsIgnoringNaN);
//                IJ.log("Size sorted: "+pixelValues.length);

                int id = (int) ((float) nPixelsIgnoringNaN*percentile);
                prctileArray[i] = pixelValues[id];
            }
            return prctileArray;
        }

    }

}