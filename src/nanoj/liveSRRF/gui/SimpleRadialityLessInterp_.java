package nanoj.liveSRRF.gui;

import ij.*;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.liveSRRF.RadialityCL;
import nanoj.liveSRRF.RadialityLessInterpCL;

public class SimpleRadialityLessInterp_ implements PlugIn {

    public String prefsHeader = this.getClass().getName();
    protected Prefs prefs = new Prefs();

    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Simple Radiality");
        gd.addNumericField("Magnification", getPrefs("magnification", 4), 0);
        gd.addNumericField("Ring-Radius", getPrefs("ringRadius", 1), 1);
        gd.addNumericField("Radiality-Sensitivity", getPrefs("sensitivity", 10), 1);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        double ringRadius = gd.getNextNumber();
        double radialitySensitivity = gd.getNextNumber();

        setPrefs("magnification", (float) magnification);
        setPrefs("ringRadius", (float) ringRadius);
        setPrefs("sensitivity", (float) radialitySensitivity);
        prefs.savePreferences();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        ImageStack imsRadiality = new ImageStack(w * magnification, h * magnification);
        ImageStack imsInterpolated = new ImageStack(w * magnification, h * magnification);
        RadialityLessInterpCL rCL = new RadialityLessInterpCL(w, h, magnification, ringRadius, radialitySensitivity);

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

    public float getPrefs(String key, float defaultValue) {
        return (float) prefs.get(prefsHeader+"."+key, defaultValue);
    }

    public void setPrefs(String key, float defaultValue) {
        prefs.set(prefsHeader+"."+key, defaultValue);
    }

}
