package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import nanoj.liveSRRF.RemoveMacroPixelArtefacts_CL;

import java.awt.*;

public class RemoveMacroPixelArtefacts_ implements PlugIn {

    String removeMPartefactVersion = "v0.2";

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        ImageStack imsAllRawData = imp.getImageStack();

        int imageWidth = imp.getWidth();
        int imageHeight = imp.getHeight();
        int imageSliceSize = imp.getNSlices();

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("liveSRRF - Macro-pixel artefact correction" + removeMPartefactVersion);
        gd.addNumericField("Magnification: ", 4, 0);
        gd.addNumericField("Number of passes: ", 1,0);
        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            return;
        }

        int magnification = (int) gd.getNextNumber();
        int nPasses = (int) gd.getNextNumber();

        IJ.log("R: "+imageWidth%magnification);
        IJ.log("R: "+imageHeight%magnification);

        if ((imageWidth%magnification != 0) || (imageHeight%magnification != 0)){
            IJ.log("Impossible magnification ! :-o");
            return;
        }

        RemoveMacroPixelArtefacts_CL artefactRemover = new RemoveMacroPixelArtefacts_CL(imageWidth, imageHeight, magnification);
        IJ.log("Number of slices to process: "+imageSliceSize);

        for (int s = 1; s <= imageSliceSize; s++) {
            IJ.log("Slice #"+s);
            artefactRemover.loadData(imsAllRawData.getProcessor(s));
            for (int i = 0; i<nPasses; i++){
                IJ.log("Pass #" + (i + 1));
                artefactRemover.correct();
            }
            artefactRemover.readBuffer();
        }

        ImageStack imsAllResults = artefactRemover.imsResults;
        ImagePlus impResults = new ImagePlus(imp.getTitle() + " - corrected", imsAllResults);
        impResults.show();
        IJ.run(impResults, "Enhance Contrast", "saturated=0.5");

        artefactRemover.release();
        IJ.log("----------------");
        IJ.log("All done.");

    }
}
