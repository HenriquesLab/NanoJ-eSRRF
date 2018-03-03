package nanoj.liveSRRF.gui;

import ij.*;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.RadialGradientConvergenceCL;

public class RadialGradientConvergence_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());

    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Radial Gradient Convergence");
        gd.addNumericField("Magnification", prefs.get("magnification", 4), 0);
        gd.addNumericField("FWHM (pixels)", prefs.get("fwhm", 3), 2);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", (float) fwhm);
        prefs.save();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        FloatProcessor fpWeightMask = new FloatProcessor(w, h);
        fpWeightMask.add(1);

        ImageStack imsPxM = new ImageStack(w * magnification, h * magnification);
        ImageStack imsRGC = new ImageStack(w * magnification, h * magnification);
        RadialGradientConvergenceCL rCL = new RadialGradientConvergenceCL(fpWeightMask, magnification, fwhm);

        for (int n=1; n<=nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fpPxMAndRGC = rCL.calculateRGC(ims.getProcessor(n), 0, 0);
            FloatProcessor fpPxM = fpPxMAndRGC[0];
            FloatProcessor fpRGC = fpPxMAndRGC[1];

            imsPxM.addSlice(fpPxM);
            imsRGC.addSlice(fpRGC);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }

        new ImagePlus(imp.getTitle()+" - Interpolated", imsPxM).show();
        new ImagePlus(imp.getTitle()+" - Radial-Gradient-Convergence", imsRGC).show();

        //rCL.showPlatforms();
        rCL.release();
    }
}
