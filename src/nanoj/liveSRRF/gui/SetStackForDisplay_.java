package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class SetStackForDisplay_ implements PlugIn {

    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();


        int nSlices = imp.getImageStack().getSize();
        int width = imp.getImageStack().getWidth();
        int height = imp.getImageStack().getHeight();
        String imageTitle = imp.getTitle();

        ImageProcessor ip;
        ImagePlus impTemp;
        String label;

        ImageStack imsInput = imp.getImageStack();
        ImageStack imsOutput = new ImageStack(width, height);

        for (int s = 1; s <= nSlices; s++) {
            ip = imsInput.getProcessor(s);
            label  = imsInput.getSliceLabel(s);

            impTemp = new ImagePlus("Slice: " + s, ip);
            IJ.run(impTemp, "Enhance Contrast", "saturated=0.5");
            IJ.run(impTemp, "8-bit", "");
//            IJ.run(impTemp, "mpl-inferno", "");

            imsOutput.addSlice(label, impTemp.getProcessor());

        }

        ImagePlus impOutput = new ImagePlus(imageTitle + " - For display", imsOutput);
//        IJ.run(impOutput, "RGB Color", "");
        impOutput.show();


    }


}
