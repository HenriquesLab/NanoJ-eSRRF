package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;

public class GetQnRmap_ implements PlugIn {

    private final boolean DEBUG = false;

    @Override
    public void run(String arg) {

        IJ.log("------------------------");
        IJ.log("------------------------");
        IJ.log("Getting QnR map...");
        String[] imageTitles = WindowManager.getImageTitles();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("eSRRF - Get QnR map ");
        gd.addChoice("RSP map: ", imageTitles, getIndexWithSubstring(imageTitles,"RSP"));
        gd.addChoice("FRC map: ", imageTitles, getIndexWithSubstring(imageTitles,"FRC"));

        gd.addCheckbox("Automatically set min and max",true);
        gd.addNumericField("FRC resolution min (in um) (default: 0.05)", 0.05, 3);
        gd.addNumericField("FRC resolution max (in um) (default: 0.2)", 0.20, 3);

        gd.addCheckbox("Remove NaNs",true);
        gd.addCheckbox("Display rescaled FRC map",false);
        gd.addCheckbox("Show logistic conversion curve",false);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        if (DEBUG) {
            IJ.log("Listing current images...");
            for (String imageTitle : imageTitles) {
                IJ.log(imageTitle);
            }
        }

        // Get the parameters from the GUI
        String RSPmapTitle = gd.getNextChoice();
        String FRCmapTitle = gd.getNextChoice();
        IJ.log("RSP map: " + RSPmapTitle);
        IJ.log("FRC map: " + FRCmapTitle);

        boolean autoResMinMax = gd.getNextBoolean();
        float Res_min = (float) gd.getNextNumber();
        float Res_max = (float) gd.getNextNumber();
        boolean removeNaNs = gd.getNextBoolean();
        boolean displayFRCrescaled = gd.getNextBoolean();
        boolean showLogisticConversionCurve = gd.getNextBoolean();

        // Get the images and prepare
        ImagePlus impRSPmap = WindowManager.getImage(RSPmapTitle);
        ImagePlus impFRCmap = WindowManager.getImage(FRCmapTitle);
        Calibration cal = impFRCmap.getCalibration();

        if (removeNaNs) {
            IJ.log("--- Removing NaNs from maps ---");
            // Run Remove NaNs 3x to make sure it's all gone.
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");

            IJ.run(impFRCmap, "Remove NaNs...", "radius=2");
            IJ.run(impFRCmap, "Remove NaNs...", "radius=2");
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");
        }

        ImageStack imsFRCmap = impFRCmap.getImageStack();
        ImageStack imsRSPmap = impRSPmap.getImageStack();

        if (autoResMinMax) {
            IJ.log("Getting min and max resolution from data...");
            float[] minmaxRes = getMinMaxfromImageStack(imsFRCmap);
            Res_min = minmaxRes[0];
            Res_max = minmaxRes[1];
        }

        // Convert FRC resolution maps via logistic conversion
        IJ.log("Logistic conversion using Min = " + Res_min + " um and Max = " + Res_max + " um...");
        ImageStack imsFRCmapConverted = logisticImageConversion(imsFRCmap, Res_min, Res_max);
        if (displayFRCrescaled) {
            ImagePlus impFRCmapConverted = new ImagePlus("FRC map Converted", imsFRCmapConverted);
            applyLUT_SQUIRREL_Errors(impFRCmapConverted);
            impFRCmapConverted.copyScale(impFRCmap);
            impFRCmapConverted.show();
            IJ.run("Maximize", "");
        }

        ImageStack imsQnRmap = calculateQnRmap(imsFRCmapConverted, imsRSPmap);
        ImagePlus impQnRmap = new ImagePlus("QnR map", imsQnRmap);
        applyLUT_SQUIRREL_Errors(impQnRmap);
        impQnRmap.copyScale(impFRCmap);
        impQnRmap.show();
        IJ.run("Maximize", "");

        if (DEBUG) {
            IJ.log("Calibration pixel width: " + cal.pixelWidth);
            IJ.log("Calibration pixel height: " + cal.pixelHeight);
            IJ.log("Calibration xOrigin: " + cal.xOrigin);
            IJ.log("Calibration yOrigin: " + cal.yOrigin);
        }

        if (showLogisticConversionCurve) displayLogisticCurve(Res_min, Res_max);

        IJ.log("-------------------");
        float[][] resultsMax = findMaximum(imsQnRmap, cal);
        for (int s = 0; s < imsQnRmap.getSize(); s++) {
            if (imsQnRmap.getSize() > 1) IJ.log("Results when using " + resultsMax[5][s] + " frames:");
            IJ.log("Max QnR value: " + resultsMax[0][s] + " at (" + (int) resultsMax[1][s] + "," + (int) resultsMax[2][s] + ")");
            IJ.log("Best Radius: " + resultsMax[3][s]);
            IJ.log("Best Sensitivity: " + resultsMax[4][s]);
            IJ.log("-------------------");
        }

        if (imsQnRmap.getSize() > 1) {
            Plot plotQnRvsFrameNumber = new Plot("QnR score vs. # frames", "# of frames", "QnR score");
            plotQnRvsFrameNumber.addPoints(resultsMax[5], resultsMax[0], Plot.LINE);
            plotQnRvsFrameNumber.show();

            Plot plotBestRadiusVsFrameNumber = new Plot("Best radius score vs. # frames", "# of frames", "Best radius");
            plotBestRadiusVsFrameNumber.addPoints(resultsMax[5], resultsMax[3], Plot.LINE);
            plotBestRadiusVsFrameNumber.show();

            Plot plotBestSensitivityVsFrameNumber = new Plot("Best sensitivity score vs. # frames", "# of frames", "Best sensitivity");
            plotBestSensitivityVsFrameNumber.addPoints(resultsMax[5], resultsMax[4], Plot.LINE);
            plotBestSensitivityVsFrameNumber.show();
        }

        IJ.run("Cascade", "");

        IJ.log("-------------------");
        IJ.log("All done.");


    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------

    public ImageStack logisticImageConversion(ImageStack ims, float min, float max){
        // This method applies a logistic conversion based on the min and max provided

        double M1 = 0.075;
        double M2 = 0.925;
        double A1 = Math.log((1-M1)/M1);
        double A2 = Math.log((1-M2)/M2);
        double x0 = (A2*(double) max - A1*(double) min)/(A2-A1);
        double k = 1/(x0-(double) max)*A1;

        ImageStack imsOut = new ImageStack();
        float[] pixels;
        float[] pixels_out;
        FloatProcessor fp;

        for (int s = 0; s < ims.getSize(); s++) {
            fp = ims.getProcessor(s+1).convertToFloatProcessor();
            pixels = (float[]) fp.getPixels();
            pixels_out = new float[pixels.length];
            for (int i = 0; i < pixels.length; i++) {
                pixels_out[i] = (float) (1 / (Math.exp(-k * (pixels[i] - x0)) + 1));
            }
            imsOut.addSlice(new FloatProcessor(ims.getWidth(), ims.getHeight(), pixels_out));
        }
        return imsOut;
    }


    public void displayLogisticCurve(float min, float max){

        double M1 = 0.075;
        double M2 = 0.925;
        double A1 = Math.log((1-M1)/M1);
        double A2 = Math.log((1-M2)/M2);
        double x0 = (A2*(double) max - A1*(double) min)/(A2-A1);
        double k = 1/(x0-(double) max)*A1;

        int n_sample = 1000;
        float[] x = new float[n_sample];
        float[] y = new float[n_sample];

        float min_display = min * (float) 0.7;
        float max_display = max * (float) 1.3;
        float delta_step = (max_display - min_display)/n_sample;

        for (int i = 0; i < n_sample; i++) {
            x[i] = min_display + delta_step*i;
            y[i] = (float) (1 / (Math.exp(-k * (x[i] - x0)) + 1));
        }

        Plot logisticCurvePlot = new Plot("Logistic curve", "Resolution", "Calibrated output (AU)");
        logisticCurvePlot.addPoints(x, y, Plot.LINE);
        logisticCurvePlot.show();

    }

    public ImageStack calculateQnRmap(ImageStack ims1, ImageStack ims2){

        // This calculates the QnR metric based on 2 maps, assuming that the images are the same dimensions.
        assert (ims1.getWidth() == ims2.getWidth() && ims1.getHeight() == ims2.getHeight() && ims1.getSize() == ims2.getSize());
        ImageStack imsOut = new ImageStack();
        FloatProcessor fp1, fp2;
        float[] pixels1, pixels2, pixels_out;

        for (int s = 0; s < ims1.getSize(); s++) {
            fp1 = ims1.getProcessor(s+1).convertToFloatProcessor();
            pixels1 = (float[]) fp1.getPixels();

            fp2 = ims2.getProcessor(s+1).convertToFloatProcessor();
            pixels2 = (float[]) fp2.getPixels();
            pixels_out = new float[pixels1.length];

            for (int i = 0; i < pixels_out.length; i++) {
                pixels_out[i] = 2 * pixels1[i] * pixels2[i] / (pixels1[i] + pixels2[i]);
            }

            imsOut.addSlice(new FloatProcessor(ims1.getWidth(), ims1.getHeight(), pixels_out));
        }
        return imsOut;
    }

    public float[][] findMaximum(ImageStack ims, Calibration cal){

        float[] pixels;
        float max_value;
        FloatProcessor fp;
        int id, x_max, y_max;

        float[][] results = new float[6][ims.getSize()];

        for (int s = 0; s < ims.getSize(); s++) {
            fp = ims.getProcessor(s+1).convertToFloatProcessor();
            pixels = (float[]) fp.getPixels();
            max_value = -Float.MAX_VALUE;
            id = 0;
            for (int i = 0; i < pixels.length; i++) {
                if (pixels[i] > max_value) {
                    max_value = pixels[i];
                    id = i;
                }
            }
            y_max = id / (ims.getWidth());
            x_max = id - y_max * ims.getWidth();

            results[0][s] = max_value;
            results[1][s] = x_max;
            results[2][s] = y_max;

            if (cal == null) {
                results[3][s] = x_max;
                results[4][s] = y_max;
                results[5][s] = s;
            }

            else {
                results[3][s] = (float) ((x_max - cal.xOrigin)*cal.pixelWidth); // Radius
                results[4][s] = (float) ((y_max - cal.yOrigin)*cal.pixelHeight); // Sensitivity
                results[5][s] = (float) ((s - cal.zOrigin)*cal.pixelDepth); // # of frames
            }

        }
        return results;

    }

    public String getIndexWithSubstring(String[] string_array, String substring){
        String entryWithSubstring = string_array[0];
        for (String s : string_array) {
            if (s.contains(substring)) entryWithSubstring = s;
        }

        return entryWithSubstring;
    }

    public static float[] getMinMaxfromImageStack(ImageStack ims){
        float[] minmax = new float[2];
        minmax[0] = Float.MAX_VALUE;
        minmax[1] = -Float.MAX_VALUE;

        double min, max;
        for (int i = 0; i < ims.getSize(); i++) {
            min = ims.getProcessor(i+1).convertToFloatProcessor().getMin();
            max = ims.getProcessor(i+1).convertToFloatProcessor().getMax();
            IJ.log("Min: "+min);
            IJ.log("Max: "+max);

            if ((float) min < minmax[0]) minmax[0] = (float) min;
            if ((float) max > minmax[1]) minmax[1] = (float) max;
        }

        return minmax;
    }


}
