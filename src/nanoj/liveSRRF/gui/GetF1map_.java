package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;
import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_FRC;


public class GetF1map_ implements PlugIn {

    private final boolean DEBUG = true;

    public void run(String arg) {

        IJ.log("------------------------");
        IJ.log("------------------------");
        IJ.log("Getting F1 map...");
        String[] imageTitles = WindowManager.getImageTitles();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("eSRRF - Get F1 map ");
        gd.addChoice("RSP map: ", imageTitles, getIndexWithSubstring(imageTitles,"RSP"));
        gd.addChoice("FRC map: ", imageTitles, getIndexWithSubstring(imageTitles,"FRC"));

        gd.addNumericField("FRC resolution min (nm) (default: 50)", 0.7, 1);
        gd.addNumericField("FRC resolution max (nm) (default: 200)", 1.3, 1);
        gd.addCheckbox("Remove NaNs",true);
        gd.showDialog();

        if (DEBUG) {
            for (int i = 0; i < imageTitles.length; i++) {
                IJ.log(imageTitles[i]);
            }
        }

        // Get the parameters from the GUI
        String RSPmapTitle = gd.getNextChoice();
        String FRCmapTitle = gd.getNextChoice();

        float Res_min = (float) gd.getNextNumber();
        float Res_max = (float) gd.getNextNumber();
        boolean removeNaNs = gd.getNextBoolean();

        // Get the images and prepare
        ImagePlus impRSPmap = WindowManager.getImage(RSPmapTitle);
        ImagePlus impFRCmap = WindowManager.getImage(FRCmapTitle);
        Calibration cal = impFRCmap.getCalibration();

        if (removeNaNs) {
            IJ.log("--- Removing NaNs from maps ---");
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");
            IJ.run(impRSPmap, "Remove NaNs...", "radius=2");

            IJ.run(impFRCmap, "Remove NaNs...", "radius=2");
            IJ.run(impFRCmap, "Remove NaNs...", "radius=2");
        }

        ImageStack imsFRCmap = impFRCmap.getImageStack();
        ImageStack imsRSPmap = impRSPmap.getImageStack();

        ImageStack imsFRCmapConverted = logisticImageConversion(imsFRCmap, Res_min, Res_max);
        ImagePlus impFRCmapConverted = new ImagePlus("FRC map Converted", imsFRCmapConverted);
        applyLUT_SQUIRREL_FRC(impFRCmapConverted);
        impFRCmapConverted.show();
        IJ.run("Maximize", "");

        ImageStack imsF1map = calculateF1map(imsFRCmapConverted, imsRSPmap);
        ImagePlus impF1map = new ImagePlus("F1 map", imsF1map);
        applyLUT_SQUIRREL_Errors(impF1map);
        impF1map.show();
        IJ.run("Maximize", "");

        float[] resultsMax = findMaximum(imsF1map);
        IJ.log("Max value: "+resultsMax[0]+ " at ("+(int) resultsMax[1]+","+(int) resultsMax[2]+")");

        if (DEBUG){
            IJ.log("Calibration pixel width: "+cal.pixelWidth);
            IJ.log("Calibration pixel height: "+cal.pixelHeight);
            IJ.log("Calibration xOrigin: "+cal.xOrigin);
            IJ.log("Calibration yOrigin: "+cal.yOrigin);
        }

        float bestRadius = (float) (((int) resultsMax[1] - cal.xOrigin)*cal.pixelWidth);
        float bestSensitivity = (float) (((int) resultsMax[2] - cal.yOrigin)*cal.pixelHeight);

        IJ.log("-------------------");
        IJ.log("Best Radius: "+bestRadius);
        IJ.log("Best Sensitivity: "+bestSensitivity);
        IJ.log("-------------------");

        IJ.run("Cascade", "");

    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------

    public ImageStack logisticImageConversion(ImageStack ims, float min, float max){
        // This method applies a logistic conversion based on the min and max provided
        FloatProcessor fp = ims.getProcessor(1).convertToFloatProcessor();
        float[] pixels = (float[]) fp.getPixels();

        double M1 = 0.075;
        double M2 = 0.925;
        double A1 = Math.log((1-M1)/M1);
        double A2 = Math.log((1-M2)/M2);
        double x0 = (A2*(double) max - A1*(double) min)/(A2-A1);
        double k = 1/(x0-(double) max)*A1;

        float[] pixels_out = new float[pixels.length];
        for (int i = 0; i < pixels.length; i++) {
            pixels_out[i] = (float) (1/(Math.exp(-k*(pixels[i]-x0))+1));
        }

        ImageStack imsOut = new ImageStack();
        imsOut.addSlice(new FloatProcessor(ims.getWidth(), ims.getHeight(), pixels_out));
        return imsOut;
    }

    public ImageStack calculateF1map(ImageStack ims1, ImageStack ims2){

        // This calculates the F1 metric based on 2 maps, assuming that the images are the same dimensions.
        assert (ims1.getWidth() == ims2.getWidth() && ims1.getHeight() == ims2.getHeight());

        FloatProcessor fp1 = ims1.getProcessor(1).convertToFloatProcessor();
        float[] pixels1 = (float[]) fp1.getPixels();

        FloatProcessor fp2 = ims2.getProcessor(1).convertToFloatProcessor();
        float[] pixels2 = (float[]) fp2.getPixels();

        float[] pixels_out = new float[pixels1.length];

        for (int i = 0; i < pixels_out.length; i++) {
            pixels_out[i] = 2*pixels1[i]*pixels2[i]/(pixels1[i] + pixels2[i]);
        }

        ImageStack imsOut = new ImageStack();
        imsOut.addSlice(new FloatProcessor(ims1.getWidth(), ims1.getHeight(), pixels_out));
        return imsOut;
    }

    public float[] findMaximum(ImageStack ims){
        FloatProcessor fp = ims.getProcessor(1).convertToFloatProcessor();
        float[] pixels = (float[]) fp.getPixels();

        float max_value = -Float.MAX_VALUE;
        int id = 0;
        for (int i = 0; i < pixels.length; i++) {
            if (pixels[i] > max_value){
                max_value = pixels[i];
                id = i;
            }
        }

        int y = id / (ims.getWidth());
        int x = id - y*ims.getWidth();

        return new float[] {max_value, x, y};

    }

    public String getIndexWithSubstring(String[] string_array, String substring){
        String entryWithSubstring = string_array[0];
        for (int i = 0; i < string_array.length; i++) {
            if (string_array[i].contains(substring)) entryWithSubstring = string_array[i];
        }

        return entryWithSubstring;
    }


}
