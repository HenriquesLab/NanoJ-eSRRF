package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.plugin.frame.RoiManager;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.array.ArrayCasting;
import nanoj.core.java.array.ArrayMath;
import nanoj.core.java.featureExtraction.ExtractRois;
import nanoj.core.java.featureExtraction.Peaks;
import nanoj.core.java.gui._BaseDialog_;
import nanoj.core.java.threading.NanoJThreadExecutor;
import nanoj.liveSRRF.GaussianFitMinimizer;

import java.awt.*;
import java.io.IOException;

import static nanoj.core.java.featureExtraction.Peaks.getROIs;
import static nanoj.core.java.featureExtraction.Peaks.populateRoiManagerWithPeaks;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;

/**
 * Created with IntelliJ IDEA.
 * User: Ricardo Henriques <paxcalpt@gmail.com>
 * Date: 02/06/17
 * Time: 14:06
 *
 * Adapted by Romain Laine, r.laine@ucl.ac.uk
 * Date: 2019-11-08
 */
public class Extract3D_PSF_ extends _BaseDialog_ {

    private RoiManager rm = null;
    private int radius, nROIs;
    private boolean showPeaks, crop, get3Dpsf;
    private float[][] peaks;
    private String fittingMethod;

    @Override
    public boolean beforeSetupDialog(String arg) {
        autoOpenImp = true;
        useSettingsObserver = true;
        return true;
    }

    @Override
    public void setupDialog() {
        gd = new NonBlockingGenericDialog("Extract PSF...");
        gd.addMessage(
                "Note: this plugin is intended to run over a image containing discrete fluorophore detections.\n");

        gd.addNumericField("Max FWHM (pixels)", getPrefs("radius", 5), 0);
        gd.addNumericField("Number of ROIs to use", getPrefs("nROIs", 100), 0);
        gd.addChoice("Fitting method: ", new String[]{"Gaussian", "Integrated gaussian"}, "Integrated gaussian");

        gd.addCheckbox("Crop PSF size to meaningful data", getPrefs("crop", true));
        gd.addCheckbox("Show detected peaks", getPrefs("showPeaks", false));
        gd.addCheckbox("Get 3D PSF", getPrefs("get3Dpsf", false));

    }

    @Override
    public boolean loadSettings() {
        // Grab data from dialog
        radius = (int) gd.getNextNumber();
        nROIs = (int) gd.getNextNumber();

        crop = gd.getNextBoolean();
        showPeaks = gd.getNextBoolean();
        get3Dpsf = gd.getNextBoolean();
        fittingMethod = gd.getNextChoice();

        setPrefs("radius", radius);
        setPrefs("nROIs", nROIs);
        setPrefs("crop", crop);

        setPrefs("showPeaks", showPeaks);
        setPrefs("get3Dpsf", get3Dpsf);

        prefs.savePreferences();

        showPreview = true;

        return true;
    }

    public void doPreview() {

        log.status("calculating preview...");

        FloatProcessor ip = imp.getProcessor().convertToFloatProcessor();
        peaks = Peaks.getPeaks(ip, nROIs, radius, 0.25f);
        rm = ExtractRois.getRoiManager();
        populateRoiManagerWithPeaks(peaks, radius, rm);
        rm.runCommand("Associate", "false");
    }

    @Override
    public void execute() throws InterruptedException, IOException {
        if (rm == null) doPreview();

        log.status("extracting ROIs...");

        ImageProcessor ip = imp.getProcessor();
        ImageStack ims = imp.getStack();
        int nSlices = imp.getNSlices();
        IJ.log("Number of slices: "+nSlices);

        Roi[] rois = getROIs(peaks, radius+2);

        // Grabbing ROIs
        ImageStack imsRoisBigger = new ImageStack(2*(radius+2)+1,2*(radius+2)+1);

        for (Roi r: rois) {
            if (r == null) continue;
            Rectangle rect = r.getBounds();
            if (rect.x+rect.width >= ip.getWidth() || rect.y+rect.height >= ip.getHeight()) {
                continue;
            }
            ip.setRoi(r);
            imsRoisBigger.addSlice(ip.crop());
        }

        int wBigger = imsRoisBigger.getWidth();
        int hBigger = imsRoisBigger.getHeight();
        int nRois = imsRoisBigger.getSize();

        double[] sigmas = new double[nRois];
        double[] amplitudes = new double[nRois];
        double[] xc = new double[nRois];
        double[] yc = new double[nRois];
        ImageStack imsRoisBiggerAligned = new ImageStack(wBigger, hBigger, nRois);

        log.status("estimating PSFs...");
        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        for (int n=0; n<nRois; n++) {
            ThreadedFitterAndRealigner t = new ThreadedFitterAndRealigner(imsRoisBigger, imsRoisBiggerAligned, n, sigmas, amplitudes, xc, yc);
            NTE.execute(t);
        }
        NTE.finish();

        int r2p1 = 2*radius+1;
        ImageStack imsRoisAligned = imsRoisBiggerAligned.crop(2,2,0,r2p1,r2p1,nRois);

        Plot plot = new Plot("Amplitude vs Sigma", "Sigma", "Amplitude");
        plot.add("cross", sigmas, amplitudes);
        plot.show();

        log.status("Calculating PSF estimate...");

        int nPixels = r2p1*r2p1;
        float[] pixelPSF = new float[nPixels];
        double amplitudesSum = ArrayMath.sum(amplitudes);

        for (int s=1; s<=nRois; s++) {
            float[] pixelsROI = (float[]) imsRoisAligned.getPixels(s);
            for (int p=0; p<nPixels;p++) pixelPSF[p] += pixelsROI[p] * amplitudes[s-1] / amplitudesSum;
            log.progress(s, nROIs);
        }
        ImageProcessor ipPSF = new FloatProcessor(r2p1, r2p1, pixelPSF);


        // iterate new model redoing weights based on correlation
        for (int i=0; i<3; i++) {
            log.status("Calculating PSF estimate... iteration "+(i+1)+"/3");

            double[] newPixelPSF = new double[nPixels];
            double weightSum = 0;

            for (int s=1; s<=nRois; s++) {
                float[] pixelsROI = (float[]) imsRoisAligned.getPixels(s);
                double weight = Math.pow(ArrayMath.calculatePPMCC(pixelsROI, pixelPSF, true), 4);
                imsRoisAligned.setSliceLabel("Weight = "+String.format("%.3g%n", weight), s);
                for (int p=0; p<nPixels;p++) newPixelPSF[p] += pixelsROI[p] * weight;
                weightSum += weight;
                log.progress(s, nROIs);
            }
            for (int p=0; p<nPixels;p++) newPixelPSF[p] /= weightSum;
            pixelPSF = ArrayCasting.doubleToFloat(newPixelPSF);
        }

        log.progress(1);
        log.status("Done...");

        if (showPeaks) {
            new ImagePlus("Detected Aligned Peaks", imsRoisAligned).show();
        }

        if (crop) {
            ImageProcessor ipPSFCropped = null;
            Rectangle roi = new Rectangle(0, 0, ipPSF.getWidth(), ipPSF.getHeight());

            float intensity = ArrayMath.getAbsAverageValue((float[]) ipPSF.getPixels()) * ipPSF.getPixelCount();
            float croppedIntensity = intensity;

            while (croppedIntensity > (intensity * 0.99f)) {
                roi.x += 1;
                roi.y += 1;
                roi.width -= 2;
                roi.height -= 2;
                ipPSF.setRoi(roi);
                ipPSFCropped = ipPSF.crop();
                croppedIntensity = ArrayMath.getAbsAverageValue((float[]) ipPSFCropped.getPixels()) * ipPSFCropped.getPixelCount();
            }

            roi.x -= 1;
            roi.y -= 1;
            roi.width += 2;
            roi.height += 2;
            ipPSF.setRoi(roi);
            ipPSF = ipPSF.crop();
        }
        ImagePlus impPSF = new ImagePlus("PSF - " + imp.getTitle(), ipPSF);
        impPSF.show();

        // Getting the 3D PSF happens here
        if (get3Dpsf){


            ImageStack imsPSF3D = new ImageStack(2*(radius+2)+1, 2*(radius+2)+1);
            for (int i = 1; i <= nSlices; i++) {
                FloatProcessor fpSlice = ims.getProcessor(i).convertToFloatProcessor();
                FloatProcessor fpAvSlice = new FloatProcessor(2*(radius+2)+1, 2*(radius+2)+1, new float[(2*(radius+2)+1)*(2*(radius+2)+1)]);

                for (int r = 0; r < nRois; r++) {

                    if (rois[r] == null) continue;
                    Rectangle rect = rois[r].getBounds();
                    if (rect.x + rect.width >= fpSlice.getWidth() || rect.y + rect.height >= fpSlice.getHeight()) {
                        continue;
                    }
                    fpSlice.setRoi(rois[r]);
                    FloatProcessor fpAlignedNormalised = interpolateFloatProcessor(fpSlice.crop().convertToFloatProcessor(), xc[r], yc[r], false); // Not normalised
                    fpAvSlice = add(fpAvSlice, fpAlignedNormalised);
                }
                imsPSF3D.addSlice(fpAvSlice);
            }

            ImagePlus imp3DPSF = new ImagePlus("3D PSF stack", imsPSF3D);
            imp3DPSF.show();
        }

    }

    class ThreadedFitterAndRealigner extends Thread {
        private final ImageStack imsOut;
        private final FloatProcessor ip;
        private final int n;
        private final double[] sigmas, amplitudes, x, y;
        private final int w, h;
        GaussianFitMinimizer fitter;

        ThreadedFitterAndRealigner(ImageStack imsIn, ImageStack imsOut, int n, double[] sigmas, double[] amplitudes, double[] x, double[] y) {
            this.imsOut = imsOut;
            this.n = n;
            this.sigmas = sigmas;
            this.x = x;
            this.y = y;
            this.amplitudes = amplitudes;
            this.ip = imsIn.getProcessor(n+1).convertToFloatProcessor();

            int model = 0; // default: Gaussian
            if (fittingMethod.equals("Gaussian")) model = GaussianFitMinimizer.GAUSSIAN;
            if (fittingMethod.equals("Integrated gaussian")) model = GaussianFitMinimizer.INTEGRATED_GAUSSIAN;

            this.fitter = new GaussianFitMinimizer(ip, 1.5, ip.getWidth()/2d, ip.getHeight()/2d, model);
            this.w = ip.getWidth();
            this.h = ip.getHeight();
        }

        @Override
        public void run() {
            fitter.calculate();
            double xc = fitter.xc;
            double yc = fitter.yc;

            double deltaX = w/2d - xc - 0.5;
            double deltaY = h/2d - yc - 0.5;

//            FloatProcessor ipAlignedNormalised = new FloatProcessor(w, h);
//
//            ip.setInterpolationMethod(ip.BICUBIC);
//
//            float vMax = -Float.MAX_VALUE;
//            float vMin = Float.MAX_VALUE;
//
//            for (int j=0; j<h; j++) {
//                for (int i=0; i<w; i++) {
//                    float v = (float) ip.getInterpolatedPixel(i-deltaX, j-deltaY); // not sure if the right direction?
//                    vMax = Math.max(vMax, v);
//                    vMin = Math.min(vMin, v);
//                    ipAlignedNormalised.setf(i, j, v);
//                }
//            }
//
//            // normalise
//            float[] pixels = (float[]) ipAlignedNormalised.getPixels();
//            for (int n=0; n<pixels.length; n++) pixels[n] = (pixels[n]-vMin)/(vMax-vMin);

            FloatProcessor ipAlignedNormalised = interpolateFloatProcessor(ip, deltaX, deltaY, true);

            this.imsOut.setProcessor(ipAlignedNormalised, this.n+1);
            this.sigmas[this.n] = fitter.sigma;
            this.amplitudes[this.n] = fitter.amplitude;
            this.x[this.n] = deltaX;
            this.y[this.n] = deltaY;
        }
    }

    public static FloatProcessor interpolateFloatProcessor(FloatProcessor fpIn, double dx, double dy, boolean normalise){

        fpIn.setInterpolationMethod(fpIn.BICUBIC);
        int w = fpIn.getWidth();
        int h = fpIn.getHeight();

        FloatProcessor fpOut = new FloatProcessor(w, h);

        float vMax = -Float.MAX_VALUE;
        float vMin = Float.MAX_VALUE;

        for (int j=0; j<h; j++) {
            for (int i=0; i<w; i++) {
                float v = (float) fpIn.getInterpolatedPixel(i-dx, j-dy); // not sure if the right direction?
                vMax = Math.max(vMax, v);
                vMin = Math.min(vMin, v);
                fpOut.setf(i, j, v);
            }
        }

        if (normalise){
            // normalise
            float[] pixels = (float[]) fpOut.getPixels();
            for (int n=0; n<pixels.length; n++) pixels[n] = (pixels[n]-vMin)/(vMax-vMin);
        }

        return fpOut;
    }

}