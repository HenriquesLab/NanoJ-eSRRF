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
    private ImageStack imsSRRF_max, imsSRRF_avg, imsSRRF_std;
    private ImagePlus impCCM = null;
    private boolean showAVG, showMAX, showSTD;

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
        gd.addNumericField("Frames per SR image", prefs.get("nFrames", 0), 0);
        gd.addChoice("Gradient estimation method", GradMethods, GradMethods[0]);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));
        gd.addCheckbox("Correct sCMOS patterning", prefs.get("correctSCMOS", false));
        gd.addMessage("-=-= Reconstructions =-=-");
        gd.addCheckbox("Show AVG reconstruction", prefs.get("showAVG", false));
        gd.addCheckbox("Show STD reconstruction", prefs.get("showSTD", true));
        gd.addCheckbox("Show MAX reconstruction", prefs.get("showMAX", false));
        gd.addCheckbox("Intensity weighting", prefs.get("intWeighting", false));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        String GradChosenMethod = gd.getNextChoice();
        boolean correctVibration = gd.getNextBoolean();
        boolean correctSCMOS = gd.getNextBoolean();

        showAVG = gd.getNextBoolean();
        showSTD = gd.getNextBoolean();
        showMAX = gd.getNextBoolean();
        boolean intWeighting = gd.getNextBoolean();


        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrames", nFrames);
        prefs.set("correctVibration", correctVibration);
        prefs.set("correctSCMOS", correctSCMOS);

        prefs.set("showAVG", showAVG);
        prefs.set("showSTD", showSTD);
        prefs.set("showMAX", showMAX);
        prefs.set("intWeighting", intWeighting);
        prefs.save();

        if (nFrames == 0) nFrames = imp.getImageStack().getSize();

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();
        int wM = ims.getWidth() * magnification;
        int hM = ims.getHeight() * magnification;
        int nPixelsM = wM * hM;

        imsSRRF_max = new ImageStack(wM, hM);
        imsSRRF_avg = new ImageStack(wM, hM);
        imsSRRF_std = new ImageStack(wM, hM);

        ImageProcessor ipRef = null; // reference slide for Cross-Correlation and vibration correction
        float[][] pixelsGRCBuffer = null; // buffer containing time-points for reconstructions

        RadialGradientConvergenceCL rCL = new RadialGradientConvergenceCL(w, h, magnification, fwhm, GradChosenMethod, intWeighting);
        ThreadedCalculateReconstructions t = null; // calculates reconstructions in parallel

        float shiftX = 0;
        float shiftY = 0;
        int counter = 0;

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
            ImageProcessor ip = imp.getProcessor();

            // Estimate vibrations
            if (correctVibration) {
                System.out.println("New reference..."+counter);
                int id = prof.startTimer();
                if (counter == 0) {
                    ipRef = ip.duplicate();
                    shiftX = 0;
                    shiftY = 0;
                }
                else {
                    float[] shift = calculateShift(ipRef, ip, 5);
                    shiftX = shift[0];
                    shiftY = shift[1];
                }
                System.out.println("Frame="+s+" shiftX="+shiftX+" shiftY="+shiftY);
                prof.recordTime("Drift Estimation", prof.endTimer(id));
            }

            // Calculate actual Radial-Gradient-Convergence
            FloatProcessor fpRGC = rCL.calculateRGC(ip, shiftX, shiftY, GradChosenMethod);
            float[] pixelsRGC = (float[]) fpRGC.getPixels();

            // Update buffer
            if (s == 1 || counter == 0) pixelsGRCBuffer = new float[min(nFrames, nSlices-s+1)][];
            pixelsGRCBuffer[counter] = pixelsRGC;

            if (counter == nFrames-1 || s == nSlices) {
                // process buffer
                if (t != null) t.finalise();
                t = new ThreadedCalculateReconstructions(pixelsGRCBuffer);
                t.start();
                counter = 0;
            }
            else counter++;
        }

        t.finalise();
        if (showMAX) new ImagePlus(imp.getTitle()+" - SRRF AVG", imsSRRF_max).show();
        if (showAVG) new ImagePlus(imp.getTitle()+" - SRRF MAX", imsSRRF_avg).show();
        if (showSTD) new ImagePlus(imp.getTitle()+" - SRRF STD", imsSRRF_std).show();

        rCL.release();
        IJ.log(prof.report()); // TODO: look into profiler since iterations add up between consective runs (by design?)
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

        //vMax = -Double.MAX_VALUE;
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

    class ThreadedCalculateReconstructions extends Thread {

        private float[][] pixelsRGCBuffer;

        public ThreadedCalculateReconstructions(float[][] pixelsRGCBuffer) {
            this.pixelsRGCBuffer = pixelsRGCBuffer;
        }

        public void run() {
            int nSlices = pixelsRGCBuffer.length;
            int nPixels = pixelsRGCBuffer[0].length;
            float[] pixelsMax = new float[nPixels];
            float[] pixelsAvg = new float[nPixels];
            float[] pixelsStd = new float[nPixels];

            for (int s=0; s<nSlices; s++) {
                for (int p=0; p<nPixels; p++) {
                    pixelsMax[p] = max(pixelsRGCBuffer[s][p], pixelsMax[p]);
                    pixelsAvg[p] += pixelsRGCBuffer[s][p] / nPixels;
                }
            }

            if (showSTD) {
                for (int s = 1; s < nSlices; s++) {
                    for (int p = 0; p < nPixels; p++) {
                        pixelsStd[p] += pow(pixelsRGCBuffer[s][p] - pixelsAvg[p], 2) / nPixels;
                    }
                }
                for (int p = 0; p < nPixels; p++) pixelsStd[p] = (float) sqrt(pixelsStd[p]);
            }

            int w = imsSRRF_max.getWidth();
            int h = imsSRRF_max.getHeight();
            if (showMAX) imsSRRF_max.addSlice(new FloatProcessor(w, h, pixelsMax));
            if (showAVG) imsSRRF_avg.addSlice(new FloatProcessor(w, h, pixelsAvg));
            if (showSTD) imsSRRF_std.addSlice(new FloatProcessor(w, h, pixelsStd));
        }

        public void finalise() {
            try {
                this.join();
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
        }
    }
}


