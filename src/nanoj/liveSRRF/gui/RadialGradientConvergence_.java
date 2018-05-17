package nanoj.liveSRRF.gui;

import ij.*;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.liveSRRF.RadialGradientConvergenceCL;

public class RadialGradientConvergence_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();


    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        // initilizaing string for gradient estimation choice
        String[] GradMethods = new String[3];
        GradMethods[0] = "3-point gradient (classic)";
        GradMethods[1] = "Robert's cross local gradient";
        GradMethods[2] = "2-point local + interpolation";

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Radial Gradient Convergence");
        gd.addNumericField("Magnification", prefs.get("magnification", 4), 0);
        gd.addNumericField("FWHM (pixels)", prefs.get("fwhm", 3), 2);
        gd.addNumericField("Sensitivity", prefs.get("sensitivity", 5), 0);
        gd.addChoice("Gradient estimation method", GradMethods, GradMethods[0]);
        gd.addCheckbox("Intensity weighting", prefs.get("intWeighting", false));
        gd.addCheckbox("Show interpolated intensity", prefs.get("showInt", false));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();
        int sensitivity = (int) gd.getNextNumber();
        String GradChosenMethod = gd.getNextChoice();
        boolean intWeighting = gd.getNextBoolean();
        boolean showInt = gd.getNextBoolean();


        IJ.log("Gradient method chosen: " + GradChosenMethod);

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("sensitivity", sensitivity);
        prefs.set("intWeighting", intWeighting);
        prefs.set("showInt", showInt);
        prefs.save();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();


        ImageStack imsRGC = new ImageStack(w * magnification, h * magnification);
        RadialGradientConvergenceCL rCL = new RadialGradientConvergenceCL(w, h, magnification, fwhm, sensitivity, GradChosenMethod, intWeighting);

        for (int n = 1; n <= nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor fpRGC = rCL.calculateRGC(ims.getProcessor(n), 0, 0, GradChosenMethod);
            imsRGC.addSlice(fpRGC);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }

        new ImagePlus(imp.getTitle() + " - Radial-Gradient-Convergence", imsRGC).show();


        if (showInt) {
            ImageStack imsInt = new ImageStack(w * magnification, h * magnification);
            for (int n = 1; n <= nSlices; n++) {
                IJ.showProgress(n, nSlices);
                FloatProcessor fpInt = rCL.calculateInt(ims.getProcessor(n), 0, 0);
                imsInt.addSlice(fpInt);

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return;
                }
            }
            new ImagePlus(imp.getTitle() + " - Interpolated intensity", imsInt).show();
        }

        prof = rCL.prof;

        //rCL.showPlatforms();
        rCL.release();
        IJ.log(prof.report()); // TODO: look into profiler since iterations add up between consective runs (by design?)

    }
}
