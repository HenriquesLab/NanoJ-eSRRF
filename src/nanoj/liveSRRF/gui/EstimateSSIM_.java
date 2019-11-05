package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import nanoj.liveSRRF.SSIMCalculator;

import java.util.Arrays;

import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;


public class EstimateSSIM_ implements PlugIn {

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

        GenericDialog gd = new GenericDialog("liveSRRF - SSIM estimator");
        gd.addChoice("Regularization method", regularizationMethods, regularizationMethods[0]);
        gd.addNumericField("Radius (in pixels)", 2, 2);
        gd.addNumericField("Frames to average for ref.", 5, 0);
        gd.addCheckbox("Temporal blocking", true);
        gd.addNumericField("Block size (frames)", 100, 0);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String regMethod = gd.getNextChoice();
        double radius = gd.getNextNumber();
        int nFramesToAverage = (int) gd.getNextNumber();
        boolean tempBlocking = gd.getNextBoolean();
        int blockSize = (int) gd.getNextNumber();

        int nBlocks = 1;
        if (tempBlocking) {
            nBlocks = (int) (float) nFrames / blockSize;
        }
        else{
            blockSize = nFrames;
        }

        for (int b = 0; b < nBlocks; b++) {

            IJ.log("Block #"+(b+1));
            ImageProcessor ipRef = ims.duplicate().getProcessor(1 + b*blockSize).convertToFloatProcessor();
            for (int i = 1; i < nFramesToAverage; i++) {
                ipRef = add(ipRef.convertToFloatProcessor(), ims.duplicate().getProcessor(i + 1 + b*blockSize).convertToFloatProcessor());
            }
            ipRef.multiply(1.0f / (double) nFramesToAverage);

            float[] regFactors = getRegularizationFactors(ims, regMethod);
            SSIMCalculator ssimCalculator = new SSIMCalculator(ipRef, radius, regFactors);
            ImageStack imsSSIM = new ImageStack(width, height, blockSize);
            for (int i = 1; i <= blockSize; i++) {
                IJ.showProgress(i, blockSize);
                ImageProcessor ip = ims.getProcessor(i + b*blockSize);
                imsSSIM.setProcessor(ssimCalculator.Calculate(ip), i);
            }

            ImagePlus impSSIM = new ImagePlus();
            impSSIM.setStack(imageTitle + " - SSIM stack (block #"+(b+1)+")", imsSSIM);
            impSSIM.show();
        }
    }

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
}