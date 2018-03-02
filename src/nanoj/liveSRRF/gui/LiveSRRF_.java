package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.RadialGradientConvergenceCL;

import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class LiveSRRF_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());

    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Radial Gradient Convergence");
        gd.addNumericField("Magnification", prefs.get("magnification", 4), 0);
        gd.addNumericField("FWHM (pixels)", prefs.get("fwhm", 3), 2);
        gd.addNumericField("Frames per SR image", prefs.get("nFrames", 0), 0);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));
        gd.addCheckbox("Correct sCMOS patterning", prefs.get("correctSCMOS", false));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        boolean correctVibration = gd.getNextBoolean();
        boolean correctSCMOS = gd.getNextBoolean();

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrames", nFrames);
        prefs.set("correctVibration", correctVibration);
        prefs.set("correctSCMOS", correctSCMOS);
        prefs.save();

        if (nFrames == 0) nFrames = imp.getImageStack().getSize();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();
        int wM = ims.getWidth() * magnification;
        int hM = ims.getHeight() * magnification;
        int nPixelsM = wM * hM;

        RadialGradientConvergenceCL rCL = new RadialGradientConvergenceCL(w, h, magnification, fwhm);

        ImageStack imsSRRF = new ImageStack(w * magnification, h * magnification);
        ImagePlus impSRRF = new ImagePlus(imp.getTitle()+" - SRRF");

        float[] pixelsSRRF_avg = null, pixelsSRRF_std = null, pixelsSRRF_max = null;

        int counter = 1;

        for (int s=1; s<=nSlices; s++) {
            //IJ.log("S="+s+" - C="+counter);
            IJ.showProgress(s, nSlices);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                rCL.release();
                return;
            }

            ImageProcessor ip = ims.getProcessor(s);
            FloatProcessor fpRGC = rCL.calculateRGC(ip, 0, 0);
            float[] pixelsRGC = (float[]) fpRGC.getPixels();

            if (s == 1 || counter == 1) {
                pixelsSRRF_max = pixelsRGC.clone();
                pixelsSRRF_avg = pixelsRGC.clone();;
                pixelsSRRF_std = new float[pixelsRGC.length];
                counter++;
            }
            else if (counter == nFrames || s == nSlices) {
                for (int p=0; p<nPixelsM; p++) pixelsSRRF_std[p] = (float) sqrt(pixelsSRRF_std[p]); // convert from variance to stddev
                imsSRRF.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_std));
                imsSRRF.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_avg));
                imsSRRF.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_max));
                pixelsSRRF_max = pixelsRGC.clone();
                pixelsSRRF_avg = pixelsRGC.clone();
                pixelsSRRF_std = new float[pixelsRGC.length];
//                counter = 2;
                counter = 1;

                // render the image
                impSRRF.setStack(imsSRRF);
                impSRRF.show();
                impSRRF.setSlice(imsSRRF.getSize());
            }
            else {
                for (int p=0; p<nPixelsM; p++) {
                    pixelsSRRF_max[p] = max(pixelsRGC[p], pixelsSRRF_max[p]);
                    pixelsSRRF_avg[p] += (pixelsRGC[p] - pixelsSRRF_avg[p]) / counter;
                    pixelsSRRF_std[p] += pow(pixelsRGC[p] - pixelsSRRF_avg[p], 2) / (counter-1);
                }
                counter++;
            }
        }
        rCL.release();

        IJ.run(impSRRF, "Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames="+(imsSRRF.getSize()/3)+" display=Composite");
        for (int s=1; s<=3; s++) {
            impSRRF.setSlice(s);
            IJ.run(impSRRF, "Enhance Contrast", "saturated=0.35"); // this does not appear to work
        }
        impSRRF.setSlice(1);
    }
}
