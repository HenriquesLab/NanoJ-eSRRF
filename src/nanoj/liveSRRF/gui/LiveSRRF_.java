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
import nanoj.core2.NanoJProfiler;
import nanoj.liveSRRF.RadialGradientConvergenceCL;

import static java.lang.Math.*;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

public class LiveSRRF_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();
    private ImageStack imsPxMAvg, imsRGCMax, imsRGCAvg, imsRGCStd;
    private ImagePlus impCCM = null;
    private boolean showAVG, showMAX, showSTD;
    private boolean correctVibration, correctSCMOS;

    private int nSlices, nPixels, nFrames, nPixelsM, w, h, wM, hM;
    private float[] pixelsPxMAvg, pixelsRGCAvg, pixelsRGCMax, pixelsRGCStd;
    private float[][] pixelsRGCBuffer;

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
        nFrames = (int) gd.getNextNumber();
        correctVibration = gd.getNextBoolean();
        correctSCMOS = gd.getNextBoolean();

        showAVG = gd.getNextBoolean();
        showSTD = gd.getNextBoolean();
        showMAX = gd.getNextBoolean();

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
        else nFrames = max(nFrames, imp.getImageStack().getSize());

        ImageStack ims = imp.getImageStack();
        nSlices = ims.getSize();
        w = ims.getWidth();
        h = ims.getHeight();
        wM = ims.getWidth() * magnification;
        hM = ims.getHeight() * magnification;
        nPixels = w * h;
        nPixelsM = wM * hM;

        imsPxMAvg = new ImageStack(wM, hM);
        imsRGCMax = new ImageStack(wM, hM);
        imsRGCAvg = new ImageStack(wM, hM);
        imsRGCStd = new ImageStack(wM, hM);

        resetArrayBuffer(); // initialised pixelsPxMAvg, pixelsRGCAvg, pixelsRGCMax, pixelsRGCStd
        float[] pixelsDarkAverage = new float[nPixels];
        float[] pixelsDarkStdDev  = new float[nPixels];
        if (correctSCMOS) {
            IJ.showStatus("Select Dark-Frames Average");
            ImagePlus impDarkAverage = IJ.openImage();
            if (impDarkAverage == null) return;
            pixelsDarkAverage = (float[]) impDarkAverage.getProcessor().convertToFloatProcessor().getPixels();
            assert (pixelsDarkAverage.length == nPixels);

            IJ.showStatus("Select Dark-Frames StdDev");
            ImagePlus impDarkStdDev = IJ.openImage();
            if (impDarkStdDev == null) return;
            pixelsDarkStdDev = (float[]) impDarkStdDev.getProcessor().convertToFloatProcessor().getPixels();
            assert (pixelsDarkStdDev.length == nPixels);
        }
        else { // if we don't have a Dark StdDev, just fill the array with 1
            for (int n=0; n<nPixels; n++) pixelsDarkStdDev[n] = 1;
        }

        FloatProcessor fpWeightMask = new FloatProcessor(w, h, pixelsDarkStdDev);
        RadialGradientConvergenceCL rCL = new RadialGradientConvergenceCL(fpWeightMask, magnification, fwhm);

        ImageProcessor fpRef = null; // reference slide for Cross-Correlation and vibration correction
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
            imp.setSlice(s);
            FloatProcessor fpFrame = imp.getProcessor().convertToFloatProcessor();
            if (correctSCMOS) {
                float[] pixels = (float[]) fpFrame.getPixels();
                for (int n=0; n<nPixels; n++) pixels[n] -= pixelsDarkAverage[n];
                //fpFrame.blurGaussian(0.5);
            }

            // Estimate vibrations
            if (correctVibration) {
                int id = prof.startTimer();
                if (counter == 0) {
                    fpRef = fpFrame.duplicate();
                    shiftX = 0;
                    shiftY = 0;
                }
                else {
                    float[] shift = calculateShift(fpRef, fpFrame, 5);
                    shiftX = shift[0];
                    shiftY = shift[1];
                }
                //System.out.println("Frame="+s+" shiftX="+shiftX+" shiftY="+shiftY);
                prof.recordTime("Drift Estimation", prof.endTimer(id));
            }

            // Calculate actual Radial-Gradient-Convergence
            FloatProcessor[] fpPxMAndRGC = rCL.calculateRGC(fpFrame, shiftX, shiftY);
            float[] pixelsPxM = (float[]) fpPxMAndRGC[0].getPixels();
            float[] pixelsRGC = (float[]) fpPxMAndRGC[1].getPixels();

            // Calculate instant projections
            for (int p=0; p<nPixelsM; p++) {
                //pixelsRGC[p] = (pixelsRGC[p] - pixelsDarkStdDev[p]) * pixelsPxM[p]; // remove offset (calculated from dark frames) from RGC
                pixelsRGC[p] = pixelsRGC[p];// * pixelsPxM[p]; // remove offset (calculated from dark frames) from RGC
                pixelsPxMAvg[p] += (pixelsPxM[p] - pixelsPxMAvg[p]) / counter;
                pixelsRGCAvg[p] += (pixelsRGC[p] - pixelsRGCAvg[p]) / counter;
                pixelsRGCMax[p] = max(pixelsRGC[p], pixelsRGCMax[p]);
            }

            // Keep a record of frames if user wants StdDev
            if (showSTD) pixelsRGCBuffer[counter-1] = pixelsRGC;

            if (counter == nFrames || s == nSlices) {
                // Re-Normalise Intensity
                float maxPxMAvg = - Float.MAX_VALUE;
                float minPxMAvg =   Float.MAX_VALUE;
                float maxRGCAvg = - Float.MAX_VALUE;
                float minRGCAvg =   Float.MAX_VALUE;
                float maxRGCMax = - Float.MAX_VALUE;
                float minRGCMax =   Float.MAX_VALUE;

                for (int p=0; p<nPixelsM; p++) {
                    maxPxMAvg = max(pixelsPxMAvg[p], maxPxMAvg);
                    minPxMAvg = min(pixelsPxMAvg[p], minPxMAvg);
                    maxRGCAvg = max(pixelsRGCAvg[p], maxRGCAvg);
                    minRGCAvg = min(pixelsRGCAvg[p], minRGCAvg);
                    maxRGCMax = max(pixelsRGCMax[p], maxRGCMax);
                    minRGCMax = min(pixelsRGCMax[p], minRGCMax);
                }

                for (int p=0; p<nPixelsM; p++) {
                    pixelsRGCAvg[p] = (pixelsRGCAvg[p] - minRGCAvg) * (maxPxMAvg - minPxMAvg) / (maxRGCAvg - minRGCAvg);
                    pixelsRGCMax[p] = (pixelsRGCMax[p] - minRGCMax) * (maxPxMAvg - minPxMAvg) / (maxRGCMax - minRGCMax);
                }

                // special case for StdDev
                if (showSTD) {
                    float maxRGCStd = - Float.MAX_VALUE;
                    float minRGCStd =   Float.MAX_VALUE;

                    for (int c=0; c<counter; c++) {
                        for (int p=0; p<nPixelsM; p++) {
                            pixelsRGCStd[p] = (float) (pow(pixelsRGCBuffer[c][p] - pixelsRGCAvg[p], 2) / counter);
                        }
                    }
                    for (int p=0; p<nPixelsM; p++) {
                        pixelsRGCStd[p] = (float) sqrt(pixelsRGCStd[p]);
                        maxRGCStd = max(pixelsRGCStd[p], maxRGCStd);
                        minRGCStd = min(pixelsRGCStd[p], minRGCStd);
                    }

                    for (int p=0; p<nPixelsM; p++)
                        pixelsRGCStd[p] = (pixelsRGCStd[p] - minRGCStd) * (maxPxMAvg - minPxMAvg) / (maxRGCStd - minRGCStd);
                }

                imsPxMAvg.addSlice(new FloatProcessor(wM, hM, pixelsPxMAvg));
                imsRGCAvg.addSlice(new FloatProcessor(wM, hM, pixelsRGCAvg));
                imsRGCMax.addSlice(new FloatProcessor(wM, hM, pixelsRGCMax));
                if (showSTD) {
                    imsRGCStd.addSlice(new FloatProcessor(wM, hM, pixelsRGCStd));
                }

                resetArrayBuffer();
            }
            else counter++;
        }

        new ImagePlus(imp.getTitle()+" - Interpolated AVG", imsPxMAvg).show();
        if (showMAX) new ImagePlus(imp.getTitle()+" - SRRF MAX", imsRGCMax).show();
        if (showAVG) new ImagePlus(imp.getTitle()+" - SRRF AVG", imsRGCAvg).show();
        if (showSTD) new ImagePlus(imp.getTitle()+" - SRRF STD", imsRGCStd).show();

        rCL.release();
        IJ.log(prof.report());
    }

    private void resetArrayBuffer() {
        pixelsPxMAvg    = new float[nPixelsM];
        pixelsRGCMax    = new float[nPixelsM];
        pixelsRGCAvg    = new float[nPixelsM];
        if (showSTD) {
            pixelsRGCStd = new float[nPixelsM];
            pixelsRGCBuffer = new float[nFrames][];
        }
    }

    private float[] calculateShift(ImageProcessor ipRef, ImageProcessor ip, int radius) {

        FloatProcessor fpCCM = (FloatProcessor) calculateCrossCorrelationMap(ipRef, ip, false);

        int windowSize = radius * 2 + 1;
        int xStart = fpCCM.getWidth() / 2 - radius;
        int yStart = fpCCM.getHeight() / 2 - radius;
        fpCCM.setRoi(xStart, yStart, windowSize, windowSize);
        fpCCM = (FloatProcessor) fpCCM.crop();

        double vMax = -Double.MAX_VALUE;
        double vMin = Double.MAX_VALUE;
        double xMax = 0;
        double yMax = 0;

        // first do coarse search for max
        for (int y = 1; y<windowSize-1; y++){
            for (int x = 1; x<windowSize-1; x++) {
                double v = fpCCM.getf(x,y);
                if (v > vMax) {
                    vMax = v;
                    xMax = x;
                    yMax = y;
                }
                vMin = min(v, vMin);
            }
        }
        //System.out.println("xMax="+xMax+" yMax="+yMax);

        // do fine search for max
        for (double y = yMax; y<yMax+1; y+=0.01){
            for (double x = xMax; x<xMax+1; x+=0.01) {
                double v = fpCCM.getBicubicInterpolatedPixel(x, y, fpCCM);
                if (v > vMax) {
                    vMax = v;
                    xMax = x;
                    yMax = y;
                }
            }
        }

        // recenter pixels
        float shiftX = (float) xMax - radius;
        float shiftY = (float) yMax - radius;

        if (impCCM == null) {
            impCCM = new ImagePlus("CCM Vibration Stabilisation", fpCCM);
            impCCM.show();
        }
        impCCM.setProcessor(fpCCM);
        impCCM.setRoi(new PointRoi(xMax+.5, yMax+.5));
        impCCM.setDisplayRange(vMin, vMax);

        return new float[] {shiftX, shiftY};
    }
}


