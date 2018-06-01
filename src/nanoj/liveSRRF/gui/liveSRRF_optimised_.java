package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import org.python.modules.math;
import nanoj.liveSRRF.liveSRRF_CL;


import java.awt.*;

import static java.lang.Math.min;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

public class liveSRRF_optimised_ implements PlugIn {

    // Basic formats
    private int magnification, nFrameForSRRF, sensitivity, frameGap, nFrameOnGPU, nSRRFframe, nSlices, width, height;
    private float fwhm, maxMemoryGPU;
    private boolean correctVibration, calculateAVG, calculateSTD, doFusion, getInterpolatedImage;
    private final int radiusCCM = 5;
    private float[] shiftX, shiftY;

    // Image formats
    private ImagePlus imp;
    private ImagePlus impCCM = null;

    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();


    @Override
    public void run(String arg) {

//        IJ.log("\\clear"); / TODO: clear the log window properly
        // Get raw data
        imp = WindowManager.getCurrentImage(); // TODO: depending on the size of data and RAM, consider using Virtual Stack load
        if (imp == null) imp = IJ.openImage();
        imp.show();

        nSlices = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("LiveSRRF");
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("FWHM (pixels, default: 2)", prefs.get("fwhm", 2), 2);
        gd.addNumericField("Sensitivity (default: 3)", prefs.get("sensitivity", 3), 0);

        gd.addNumericField("# frames for SRRF (0 = auto)", prefs.get("nFrameForSRRF", 0), 0);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("calculateSTD", false));
        gd.addCheckbox("Fusion reconstruction (default: off)", prefs.get("doFusion", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("getInterpolatedImage", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addNumericField("Gap between SR frame (frames, default: 50)", prefs.get("frameGap", 50), 0);
        gd.addMessage("Rolling analysis may lead to long computation times.");

        gd.addMessage("-=-= Memory =-=-\n", headerFont);
        gd.addNumericField("Maximum amount of memory on GPU (MB, default: 1000)", prefs.get("maxMemoryGPU", 500), 2);
        gd.addMessage("Giving SRRF access to a lot of memory speeds up the reconstruction\n" +
                "but may slow down the graphics card for your Minecraft game that you have \n" +
                "running in parallel.");


        MyDialogListener dl = new MyDialogListener(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        // Save last user entries
        savePreferences();

        shiftX = new float[nFrameForSRRF];
        shiftY = new float[nFrameForSRRF];

        if (calculateAVG || doFusion) {
            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG projection)");
            impSRRFavg.copyScale(imp); // make sure we copy the pixel sizes correctly across
            Calibration cal = impSRRFavg.getCalibration();
            cal.pixelWidth /= magnification;
            cal.pixelHeight /= magnification;

            ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);

        }

        if (calculateSTD || doFusion) {
            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD projection)");
            impSRRFstd.copyScale(imp); // make sure we copy the pixel sizes correctly across
            Calibration cal = impSRRFstd.getCalibration();
            cal.pixelWidth /= magnification;
            cal.pixelHeight /= magnification;

            ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);
        }

        if (doFusion) {
            ImagePlus impSRRFfusion = new ImagePlus(imp.getTitle() + " - liveSRRF (fusion)");
            impSRRFfusion.copyScale(imp); // make sure we copy the pixel sizes correctly across
            Calibration cal = impSRRFfusion.getCalibration();
            cal.pixelWidth /= magnification;
            cal.pixelHeight /= magnification;

            ImageStack imsSRRFfusion = new ImageStack(width * magnification, height * magnification);
        }


        IJ.log("-------------------------------------");
        IJ.log("liveSRRF - SARF");
        IJ.log("Magnification: " + magnification);
        IJ.log("# frames for SRRF: " + nFrameForSRRF);
        IJ.log("# frames gap: " + frameGap);
        IJ.log("# frames on GPU: " + nFrameOnGPU);
        IJ.log("GPU memory usage: " + Math.round(predictMemoryUsed(nFrameOnGPU)[0]) + " MB");
        IJ.log("RAM usage: " + Math.round(predictMemoryUsed(nFrameOnGPU)[1]) + " MB");
        IJ.log("# reconstructed frames: " + nSRRFframe);

        int indexStartSRRFframe;
        int indexStartGPUloadframe;
        int nFrameToLoad = nFrameOnGPU;
        int counterGPUframe = 0;
        int counterRawDataLoad = 0;
        liveSRRF_CL liveSRRF = new liveSRRF_CL(width, height, magnification, fwhm, sensitivity, nFrameOnGPU, nFrameForSRRF);
        liveSRRF.loadRawDataGPUbuffer(imp, 0, nFrameOnGPU);

        for (int r = 1; r <= nSRRFframe; r++) {
            IJ.log("--------");
            IJ.log("SRRF frame: " + r);
            indexStartSRRFframe = (r - 1) * frameGap + 1;
            IJ.log("Stack index start: " + indexStartSRRFframe);
            if (correctVibration) calculateShiftArray(indexStartSRRFframe);

            for (int s = 1; s <= nFrameForSRRF; s++) {
                // calculateSRRF image
                counterGPUframe++;
                if (counterGPUframe > nFrameOnGPU) {
                    counterGPUframe = 0;
                    counterRawDataLoad++;
                    nFrameToLoad = min(nFrameOnGPU, nSlices - nFrameOnGPU * counterRawDataLoad); //TODO: exclude the last few frames that are not analysed
                    liveSRRF.loadRawDataGPUbuffer(imp, nFrameOnGPU * counterRawDataLoad, nFrameToLoad);
                }

            }


        }

        liveSRRF.release(); // Release the GPU!!!
        IJ.log("-------------------------------------");
        IJ.log("Thank you for your custom on this beautiful day !");


    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------

    //    --- Grab settings ---
    private boolean grabSettings(GenericDialog gd) {

        boolean goodToGo = false;


        magnification = (int) gd.getNextNumber();
        fwhm = (float) gd.getNextNumber();
        sensitivity = (int) gd.getNextNumber();
        nFrameForSRRF = (int) gd.getNextNumber();

        correctVibration = gd.getNextBoolean();

        calculateAVG = gd.getNextBoolean();
        calculateSTD = gd.getNextBoolean();
        doFusion = gd.getNextBoolean();
        getInterpolatedImage = gd.getNextBoolean();

        frameGap = (int) gd.getNextNumber();
        maxMemoryGPU = (int) gd.getNextNumber();

        if (nFrameForSRRF == 0) nFrameForSRRF = nSlices;
        nFrameForSRRF = min(nSlices, nFrameForSRRF);
        nSRRFframe = (int) ((float) (nSlices - nFrameForSRRF) / (float) frameGap) + 1;

        if (frameGap == 0) frameGap = nFrameForSRRF;

        int maxnFrameOnGPU = 0;
        double[] memUsed = new double[2];
        while (memUsed[0] < (double) maxMemoryGPU) {
            maxnFrameOnGPU = maxnFrameOnGPU + 1;
            memUsed = predictMemoryUsed(maxnFrameOnGPU);
        }

        maxnFrameOnGPU = maxnFrameOnGPU - 1;
        if (maxnFrameOnGPU > 0) {
            goodToGo = true;
            nFrameOnGPU = (int) math.ceil((float) nSlices / (float) maxnFrameOnGPU);
            nFrameOnGPU = (int) math.ceil(((float) nSlices / (float) nFrameOnGPU));
            IJ.showStatus("liveSRRF - Number of frames on GPU: " + (nFrameOnGPU));
        } else {
            memUsed = predictMemoryUsed(1);
            IJ.showStatus("liveSRRF - Minimum GPU memory: " + Math.round(memUsed[0]) + "MB");
        }


        return goodToGo;

    }

    //    --- Dialog listener ---
    class MyDialogListener implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {

            boolean goodToGo = grabSettings(gd);
            return goodToGo;

        }
    }

    // --- Predict memory usage in MB ---
    private double[] predictMemoryUsed(int nFrameOnGPU) {

        double[] memUsed = new double[2];

        // Memory on GPU ----
        memUsed[0] = 0;
        memUsed[0] += width * height * nFrameOnGPU; // clBufferPx
        memUsed[0] += nFrameOnGPU; // clBufferShiftX
        memUsed[0] += nFrameOnGPU; // clBufferShiftY
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGx
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGy
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGxInt
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGyInt
        memUsed[0] += 2 * width * height * magnification * magnification; // clBufferRGC
        if (getInterpolatedImage) memUsed[0] += width * height * magnification * magnification; // clBufferInt

        // Memory on CPU ----
        memUsed[1] = 0;
        memUsed[1] += width * height * nSlices; // clBufferPx
        memUsed[1] += 2 * nSRRFframe * width * height * magnification * magnification; // clBufferSRRF // TODO: estimating RAM requirements properly
        if (getInterpolatedImage) memUsed[1] += width * height * magnification * magnification; // clBufferInt

        memUsed[0] *= Float.SIZE / 8000000d;
        memUsed[1] *= Float.SIZE / 8000000d;

        return memUsed;
    }


    // --- Calculate shift using Cross-correlation matrix ---
    private float[] calculateShift(ImageProcessor ipRef, ImageProcessor ipData) {

        FloatProcessor fpCCM = (FloatProcessor) calculateCrossCorrelationMap(ipRef, ipData, false);

        int windowSize = radiusCCM * 2 + 1;
        int xStart = fpCCM.getWidth() / 2 - radiusCCM;
        int yStart = fpCCM.getHeight() / 2 - radiusCCM;
        fpCCM.setRoi(xStart, yStart, windowSize, windowSize);
        fpCCM = (FloatProcessor) fpCCM.crop();

        double vMax = -Double.MAX_VALUE;
        double vMin = Double.MAX_VALUE;
        double xMax = 0;
        double yMax = 0;

        // first do coarse search for max
        for (int y = 1; y < windowSize - 1; y++) {
            for (int x = 1; x < windowSize - 1; x++) {
                double v = fpCCM.getf(x, y);
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
        for (double y = yMax; y < yMax + 1; y += 0.01) {
            for (double x = xMax; x < xMax + 1; x += 0.01) {
                double v = fpCCM.getBicubicInterpolatedPixel(x, y, fpCCM);
                if (v > vMax) {
                    vMax = v;
                    xMax = x;
                    yMax = y;
                }
            }
        }

        // recenter pixels
        float shiftX = (float) xMax - radiusCCM;
        float shiftY = (float) yMax - radiusCCM;

        if (impCCM == null) {
            impCCM = new ImagePlus("CCM Vibration Stabilisation", fpCCM);
            impCCM.show();
        }

        impCCM.setProcessor(fpCCM);
        impCCM.setRoi(new PointRoi(xMax + .5, yMax + .5));
        impCCM.setDisplayRange(vMin, vMax);

        return new float[]{shiftX, shiftY};
    }


    // -- Calculate shift Array using Cross-correlation matrix --
    private void calculateShiftArray(int indexStart) {

        imp.setSlice(indexStart);
        ImageProcessor ipRef = imp.getProcessor();
        ImageProcessor ipData;

        for (int s = 0; s < nFrameForSRRF; s++) {

            // Check if user is cancelling calculation
            IJ.showProgress(s, nFrameForSRRF);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }

            // Grab the new frame from the list
            imp.setSlice(s + indexStart);
            ipData = imp.getProcessor();

            // Estimate vibrations
            int id = prof.startTimer();
            float[] shift = calculateShift(ipRef, ipData);
            shiftX[s] = shift[0];
            shiftY[s] = shift[1];

            System.out.println("Frame=" + s + " shiftX=" + shiftX[s] + " shiftY=" + shiftY[s]);
            prof.recordTime("Drift Estimation", prof.endTimer(id));


        }

    }

    // -- Save user's last parameter set --
    private void savePreferences() {

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrameForSRRF", nFrameForSRRF);
        prefs.set("sensitivity", sensitivity);
        prefs.set("correctVibration", correctVibration);

        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateSTD", calculateSTD);
        prefs.set("doFusion", doFusion);
        prefs.set("getInterpolatedImage", getInterpolatedImage);

        prefs.set("frameGap", frameGap);
        prefs.set("maxMemoryGPU", maxMemoryGPU);

        prefs.save();

    }


}
