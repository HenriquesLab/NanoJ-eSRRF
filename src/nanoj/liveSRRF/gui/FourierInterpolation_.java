package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import static nanoj.liveSRRF.FourierInterpolation.*;

public class FourierInterpolation_ implements PlugIn { // TODO: check what happens when using any sized image
    public void run(String arg) {
        IJ.log("\\Clear");  // Clear the log window

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        String[] transformList = new String[2];
        int n = 0;
        transformList[n++] = "Interpolation";
        transformList[n++] = "Mirror padding";

        GenericDialog gd = new GenericDialog("FHT tools");
        gd.addChoice("Transform: ", transformList, transformList[1]);
        gd.addNumericField("Interpolation factor", 4, 0);
        gd.addCheckbox("Do mirror padding?", true);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String chosenTransform = gd.getNextChoice();
        int interpolationFactor = (int) gd.getNextNumber();
        boolean doMirrorPadding = gd.getNextBoolean();
        ImageStack ims = imp.getStack();
        int nFrames = ims.getSize();

        ImageStack imsInt = new ImageStack();
        ImageStack imsPadded = new ImageStack();


        if (chosenTransform.equals("Interpolation")) {
            for (int i = 1; i <= nFrames; i++) {
//            // Mirror padding of the original image
                FloatProcessor fp = ims.getProcessor(i).duplicate().convertToFloatProcessor();
                FloatProcessor fpInt = fhtSpaceInterpolation(fp, interpolationFactor, doMirrorPadding);

                if (doMirrorPadding) {
                    FloatProcessor fpPadded = mirrorPaddingEvenSquare(fp);
                    imsPadded.addSlice(fpPadded);
                }

                imsInt.addSlice(fpInt);
            }

            ImagePlus impInt = new ImagePlus("FHT-interpolated image", imsInt);
            impInt.show();
            IJ.run(impInt, "Enhance Contrast", "saturated=0.35");

            if (doMirrorPadding) {
                ImagePlus impPadded = new ImagePlus("Mirror-padded image", imsPadded);
                impPadded.show();
                IJ.run(impPadded, "Enhance Contrast", "saturated=0.35");
            }
        }

        if (chosenTransform.equals("Mirror padding")) {
            for (int i = 1; i <= nFrames; i++) {
                FloatProcessor fp = ims.getProcessor(i).duplicate().convertToFloatProcessor();
                FloatProcessor fpPadded = mirrorPaddingEvenSquare(fp);
                imsPadded.addSlice(fpPadded);
            }

            ImagePlus impPadded = new ImagePlus("Padded image", imsPadded);
            impPadded.show();
            IJ.run(impPadded, "Enhance Contrast", "saturated=0.35");

        }


    }
}