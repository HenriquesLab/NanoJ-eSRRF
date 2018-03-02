package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.RadialGradientConvergenceCL;

import static java.lang.Math.*;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

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
        gd.addMessage("-=-= Reconstructions =-=-");
        gd.addCheckbox("Show AVG reconstruction", prefs.get("showAVG", false));
        gd.addCheckbox("Show STD reconstruction", prefs.get("showSTD", true));
        gd.addCheckbox("Show MAX reconstruction", prefs.get("showMAX", false));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        boolean correctVibration = gd.getNextBoolean();
        boolean correctSCMOS = gd.getNextBoolean();

        boolean showAVG = gd.getNextBoolean();
        boolean showSTD = gd.getNextBoolean();
        boolean showMAX = gd.getNextBoolean();

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrames", nFrames);
        prefs.set("correctVibration", correctVibration);
        prefs.set("correctSCMOS", correctSCMOS);

        prefs.set("showAVG", showAVG);
        prefs.set("showSTD", showSTD);
        prefs.set("showMAX", showMAX);
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

        ImageStack imsSRRF_max = new ImageStack(w * magnification, h * magnification);
        ImageStack imsSRRF_avg = new ImageStack(w * magnification, h * magnification);
        ImageStack imsSRRF_std = new ImageStack(w * magnification, h * magnification);
        ImageStack imsCCM     = new ImageStack(11, 11);
        ImagePlus impSRRF_max = new ImagePlus(imp.getTitle()+" - SRRF AVG");
        ImagePlus impSRRF_avg = new ImagePlus(imp.getTitle()+" - SRRF MAX");
        ImagePlus impSRRF_std = new ImagePlus(imp.getTitle()+" - SRRF STD");
        ImagePlus impCCM      = new ImagePlus(imp.getTitle()+" - CCM");

        ImageProcessor ipRef = null; // reference slide for Cross-Correlation and vibration correction

        float[] pixelsSRRF_avg = null, pixelsSRRF_std = null, pixelsSRRF_max = null;

        float shiftX = 0;
        float shiftY = 0;
        int counter = 1;

        for (int s=1; s<=nSlices; s++) {
            // Check if user is cancelling calculation
            IJ.showProgress(s, nSlices);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                rCL.release();
                return;
            }

            // Grab the new frame from the list
            ImageProcessor ip = ims.getProcessor(s);

            // Estimate vibrations
            if (correctVibration) {
                if (counter == 1) ipRef = ip;

                FloatProcessor fpCCM = (FloatProcessor) calculateCrossCorrelationMap(ipRef, ip, false);
                // lets assume no more than 5 pixels shift
                int xStart = fpCCM.getWidth() / 2 - 5;
                int yStart = fpCCM.getHeight() / 2 - 5;
                fpCCM.setRoi(xStart, yStart, 11, 11);
                fpCCM = (FloatProcessor) fpCCM.crop();
                double vMax = -Double.MAX_VALUE;
                double vMin = Double.MAX_VALUE;
                for (double y = 1; y<10; y+=0.1){
                    for (double x = 1; x<10; x+=0.1) {
                        double v = fpCCM.getBicubicInterpolatedPixel(x, y, fpCCM);
                        if (v > vMax) {
                            vMax = v;
                            shiftX = (float) x - 5;
                            shiftY = (float) y - 5;
                        }
                        vMin = min(v, vMin);
                    }
                }
                System.out.println("Frame="+s+" shiftX="+shiftX+" shiftY="+shiftY);
                PointRoi r = new PointRoi(shiftX+5.5, shiftX+5.5);
                imsCCM.addSlice(fpCCM);
                impCCM.setStack(imsCCM);
                impCCM.show();
                impCCM.setSlice(imsCCM.getSize());
                impCCM.setRoi(r);
                impCCM.setDisplayRange(vMin, vMax);
            }

            // Calculate actual Radial-Gradient-Convergence
            FloatProcessor fpRGC = rCL.calculateRGC(ip, shiftX, shiftY);
            float[] pixelsRGC = (float[]) fpRGC.getPixels();

            if (s == 1 || counter == 1) {
                pixelsSRRF_max = pixelsRGC.clone();
                pixelsSRRF_avg = pixelsRGC.clone();;
                pixelsSRRF_std = new float[pixelsRGC.length];
                counter++;
            }
            else if (counter == nFrames || s == nSlices) {
                for (int p=0; p<nPixelsM; p++) pixelsSRRF_std[p] = (float) sqrt(pixelsSRRF_std[p]); // convert from variance to stddev
                imsSRRF_std.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_std));
                imsSRRF_avg.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_avg));
                imsSRRF_max.addSlice(new FloatProcessor(wM, hM, pixelsSRRF_max));
                pixelsSRRF_max = pixelsRGC.clone();
                pixelsSRRF_avg = pixelsRGC.clone();
                pixelsSRRF_std = new float[pixelsRGC.length];
                counter = 1;

                // render the image
                if (showAVG) {
                    impSRRF_avg.setStack(imsSRRF_avg);
                    impSRRF_avg.show();
                    impSRRF_avg.setSlice(imsSRRF_avg.getSize());
                }
                if (showSTD) {
                    impSRRF_std.setStack(imsSRRF_std);
                    impSRRF_std.show();
                    impSRRF_std.setSlice(imsSRRF_std.getSize());
                }
                if (showMAX) {
                    impSRRF_max.setStack(imsSRRF_max);
                    impSRRF_max.show();
                    impSRRF_max.setSlice(imsSRRF_max.getSize());
                }
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

//        IJ.run(impSRRF, "Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=1 frames="+(imsSRRF.getSize()/3)+" display=Composite");
//        for (int s=1; s<=3; s++) {
//            impSRRF.setSlice(s);
//            IJ.run(impSRRF, "Enhance Contrast", "saturated=0.35"); // this does not appear to work
//        }
//        impSRRF.setSlice(1);
    }
}
