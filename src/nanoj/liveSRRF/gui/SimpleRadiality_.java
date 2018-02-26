package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.liveSRRF.RadialityCL;

public class SimpleRadiality_ implements PlugIn {

    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Simple Radiality");
        gd.addNumericField("Magnification", 4, 0);
        gd.addNumericField("Ring-Radius", 1, 1);
        gd.addNumericField("Radiality-Sensitivity", 10, 1);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        double ringRadius = gd.getNextNumber();
        double radialitySensitivity = gd.getNextNumber();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        ImageStack imsRadiality = new ImageStack(w * magnification, h * magnification);
        ImageStack imsInterpolated = new ImageStack(w * magnification, h * magnification);
        RadialityCL rCL = new RadialityCL(w, h, magnification, ringRadius, radialitySensitivity);

        for (int n=1; n<=nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fps = rCL.calculateRadiality(ims.getProcessor(n), 0, 0);
            FloatProcessor fpRad = fps[0];
            FloatProcessor fpInterpolated = fps[1];

            imsRadiality.addSlice(fpRad);
            imsInterpolated.addSlice(fpInterpolated);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }

        new ImagePlus(imp.getTitle()+" - Radiality", imsRadiality).show();
        new ImagePlus(imp.getTitle()+" - Interpolated", imsInterpolated).show();

        //rCL.showPlatforms();
        rCL.release();
    }
}
