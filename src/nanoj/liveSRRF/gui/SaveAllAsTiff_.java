package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.PlugIn;

public class SaveAllAsTiff_ implements PlugIn {

    @Override
    public void run(String s) {

        IJ.log("------------------------");
        String[] imageTitles = WindowManager.getImageTitles();

        String pathToDisk = IJ.getDirectory("Choose where to save the images...");
        IJ.log("Saving all open images as Tiff...");

        for (String imageTitle : imageTitles) {
            ImagePlus imp = WindowManager.getImage(imageTitle);
            IJ.saveAs(imp, "Tiff", pathToDisk + imageTitle);
        }

        IJ.log("All done.");
    }
}
