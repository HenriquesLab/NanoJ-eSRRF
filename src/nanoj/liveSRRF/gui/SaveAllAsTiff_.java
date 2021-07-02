package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;

public class SaveAllAsTiff_ implements PlugIn {

    @Override
    public void run(String s) {

        IJ.log("------------------------");
        IJ.log("------------------------");
        IJ.log("Getting F1 map...");
        String[] imageTitles = WindowManager.getImageTitles();

        String pathToDisk = IJ.getDirectory("Choose where to save the images...");

        for (int i = 0; i < imageTitles.length; i++) {
            ImagePlus imp = WindowManager.getImage(imageTitles[i]);
            IJ.saveAs(imp, "Tiff", pathToDisk+imageTitles[i]);
        }

        IJ.log("All done.");
    }
}
