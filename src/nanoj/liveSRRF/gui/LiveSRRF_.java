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
import nanoj.core2.NanoJUsageTracker;
import nanoj.liveSRRF.liveSRRF_CL;

import java.awt.*;

import static java.lang.Math.*;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;
import static nanoj.liveSRRF.liveSRRF_CL.*;

public class LiveSRRF_ implements PlugIn {

    private Font headerFont = new Font("Arial", Font.BOLD, 16);

    private String user = "FijiUpdater";
    private String version = "20180418-" + user;

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJUsageTracker tracker = new NanoJUsageTracker("NanoJ-LiveSRRF", version, "UA-61590656-4");
    private NanoJProfiler prof = new NanoJProfiler();

    private ImageStack imsSRRF;
    private ImagePlus impCCM = null;

    private ImagePlus imp;
    private int magnification, nFrameSRRF, sensitivity, frameGap, nFrameOnGPU, nBufferFrames;
    private float fwhm, maxMemoryGPU;
    private boolean correctVibration, showAVGReconstruction, showSTDReconstruction, showFusionReconstruction;

    @Override
    public void run(String arg) {
        tracker.logUsage("LiveSRRF_");
        nFrameOnGPU = 1; // initial start

        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("LiveSRRF - Unfinished...");
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("FWHM (pixels, default: 2)", prefs.get("fwhm", 2), 2);
        gd.addNumericField("Sensitivity (default: 3)", prefs.get("fwhm", 3), 0);

        gd.addNumericField("Frames_per_time-point (0 - auto)", prefs.get("nFrameSRRF", 0), 0);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("Show AVG reconstruction (default: on)", prefs.get("showAVGReconstruction", true));
        gd.addCheckbox("Show STD reconstruction (default: off)", prefs.get("showSTDReconstruction", false));
        gd.addCheckbox("Show Fusion reconstruction (default: off)", prefs.get("showFusionReconstruction", false));
        gd.addCheckbox("Show Wide-field interpolation (default: off)", prefs.get("showWFInterpolation", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addNumericField("Gap between SR frame (frames, default: 50)", prefs.get("frameGap_RA", 50), 0);

        gd.addMessage("-=-= Memory =-=-\n", headerFont);
        gd.addNumericField("Maximum amount of memory on GPU (MB, default: 1000)", prefs.get("maxMemGPU", 1000), 3);

        gd.addMessage("-=-= Advice =-=-\n", headerFont);
        gd.addMessage(
                "liveSRRF is a GPU resources-hog, if it fails to run consider\n" +
                        "reducing the magnification. For larger temporal datasets,\n" +
                        "consider doing batches of 100 frames at a time, then do a \n" +
                        "a mean projecting of them.");

        gd.addMessage("-=-= CONFIDENTIAL =-=-\n", headerFont);
        gd.addMessage(
                "This is a preview of the 'in development' version of the \n" +
                        "liveSRRF engine developed by the Henriques lab @ UCL.\n" +
                        "It is only meant to be used by researchers who received\n" +
                        "a direct email from Ricardo Henriques.");

        MyDialogListener dl = new MyDialogListener(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        ImageStack ims = imp.getImageStack();

        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        int wM = w * magnification;
        int hM = h * magnification;

        ImagePlus impSRRF = new ImagePlus(imp.getTitle() + " - liveSRRF");
        impSRRF.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impSRRF.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;

        imsSRRF = new ImageStack(wM, hM);

        ImageStack imsRawDataBuffer = new ImageStack(w, h);
        ImageProcessor ip;

        boolean firstTime = true;

        ImageProcessor ipRef = null; // reference slide for Cross-Correlation and vibration correction
        liveSRRF_CL liveSRRF = new liveSRRF_CL(w, h, nFrameOnGPU, magnification, fwhm, sensitivity);

        float[] shiftX = new float[nFrameSRRF];
        float[] shiftY = new float[nFrameSRRF];

        int nFrameOnGPUcounter = 0;

//        //////////////////////////////////////
//        // !!! MAIN LOOP THROUGH FRAMES !!! //
//        //////////////////////////////////////
//
        for (int s = 1; s <= nSlices; s++) {
            // Check if user is cancelling calculation
            IJ.showProgress(s, nSlices);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                liveSRRF.release();
                return;
            }

            // Grab the new frame from the list
            imp.setSlice(s);
            ip = imp.getProcessor();
            imsRawDataBuffer.addSlice(ip);

            // Estimate vibrations
            if (correctVibration) {
                System.out.println("New reference..." + nFrameOnGPUcounter);
                int id = prof.startTimer();
                if (nFrameOnGPUcounter == 0) {
                    ipRef = ip.duplicate();
                    shiftX[nFrameOnGPUcounter] = 0;
                    shiftY[nFrameOnGPUcounter] = 0;
                } else {
                    float[] shift = calculateShift(ipRef, ip, 5);
                    shiftX[nFrameOnGPUcounter] = shift[0];
                    shiftY[nFrameOnGPUcounter] = shift[1];
                }
                System.out.println("Frame=" + s + " shiftX=" + shiftX[nFrameOnGPUcounter] + " shiftY=" + shiftY[nFrameOnGPUcounter]);
                prof.recordTime("Drift Estimation", prof.endTimer(id));
            }

//            if (counter == nFrameSRRF - 1 || s == nSlices) {
//                int id = prof.startTimer();
//                ImageStack imsResults = liveSRRF.calculateSRRF(imsRawDataBuffer, shiftX, shiftY);
//                nReconstructions = imsResults.getSize();
//                if (!showAllReconstructions) {
//                    imsSRRF.addSlice(imsResults.getProcessor(nReconstructions));
//                    imsSRRF.setSliceLabel(liveSRRF.reconstructionLabel.get(nReconstructions - 1), imsSRRF.getSize());
//                } else {
//                    for (int r = 1; r <= imsResults.getSize(); r++) {
//                        imsSRRF.addSlice(imsResults.getProcessor(r));
//                        imsSRRF.setSliceLabel(liveSRRF.reconstructionLabel.get(r - 1), imsSRRF.getSize());
//                    }
//                }
//                if (firstTime) {
//                    impSRRF.setStack(imsSRRF);
//                    impSRRF.show();
//                    impSRRF.setSlice(imsSRRF.getSize());
//                    IJ.run(impSRRF, "Enhance Contrast", "saturated=0.35");
//                    firstTime = false;
//                } else {
//                    impSRRF.setSlice(imsSRRF.getSize());
//                }
//
//                // reset buffers
//                imsRawDataBuffer = new ImageStack(w, h);
//                counter = 0;
//                prof.recordTime("full SRRF-frame calculation", prof.endTimer(id));
//            } else counter++;
        }

        liveSRRF.release(); // Release the GPU!!!
        IJ.log(prof.report());

        // Show final rendering...
        impSRRF.setStack(imsSRRF);
        IJ.run(impSRRF, "Enhance Contrast", "saturated=0.5");
        impSRRF.setTitle(imp.getTitle() + " - liveSRRF");


    }


//    -------------------------------------------------------------------------------------
//    -------------------------- Here lies dragons and functions --------------------------
//    -------------------------------------------------------------------------------------

    //    --- Grab settings ---
    private boolean grabSettings(GenericDialog gd) {

        boolean goodToGo = false;
        magnification = (int) gd.getNextNumber();
        fwhm = (float) gd.getNextNumber();
        nFrameSRRF = (int) gd.getNextNumber();
        sensitivity = (int) gd.getNextNumber();
        correctVibration = gd.getNextBoolean();

        showAVGReconstruction = gd.getNextBoolean();
        showSTDReconstruction = gd.getNextBoolean();
        showFusionReconstruction = gd.getNextBoolean();

        frameGap = (int) gd.getNextNumber();
        maxMemoryGPU = (int) gd.getNextNumber();

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrameSRRF", nFrameSRRF);
        prefs.set("sensitivity", sensitivity);
        prefs.set("correctVibration", correctVibration);

        prefs.set("showAVGReconstruction", showAVGReconstruction);
        prefs.set("showSTDReconstruction", showSTDReconstruction);
        prefs.set("showFusionReconstruction", showFusionReconstruction);
        prefs.set("showWFInterpolation", showFusionReconstruction);

        prefs.set("frameGap", frameGap);
        prefs.set("maxMemoryGPU", maxMemoryGPU);

        prefs.save();

        nBufferFrames = (int) ((float) nFrameSRRF / (float) frameGap);

        if (nFrameSRRF == 0) nFrameSRRF = imp.getImageStack().getSize();
        nFrameSRRF = min(imp.getImageStack().getSize(), nFrameSRRF);

        nFrameOnGPU = 0;
        double[] memUsed = new double[2];
        while (memUsed[0] < (double) maxMemoryGPU) {
            nFrameOnGPU = nFrameOnGPU + frameGap;
            memUsed = predictMemoryUsed(nFrameOnGPU);
            IJ.showStatus("liveSRRF - Number of frames on GPU: " + (nFrameOnGPU - frameGap));
        }

        nFrameOnGPU = nFrameOnGPU - frameGap;
        if (nFrameOnGPU > 0) goodToGo = true;
        else {
            memUsed = predictMemoryUsed(frameGap);
            IJ.showStatus("liveSRRF - Minimum GPU memory: " + Math.round(memUsed[1]) + "MB");
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

    // --- Calculate shift using Cross-correlation matrix ---
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
        float shiftX = (float) xMax - radius;
        float shiftY = (float) yMax - radius;

        if (impCCM == null) {
            impCCM = new ImagePlus("CCM Vibration Stabilisation", fpCCM);
            impCCM.show();
        }
        impCCM.setProcessor(fpCCM);
        impCCM.setRoi(new PointRoi(xMax + .5, yMax + .5));
        impCCM.setDisplayRange(vMin, vMax);

        return new float[]{shiftX, shiftY};
    }

    // --- Predict memory usage in MB ---
    private double[] predictMemoryUsed(int nFrameOnGPU) {

        double[] memUsed = new double[2];

        float width = imp.getImageStack().getWidth();
        float height = imp.getImageStack().getHeight();
        float nSlices = imp.getImageStack().getSize();

        // Memory on GPU ----
        memUsed[0] = 0;
        memUsed[0] += width * height * nFrameOnGPU; // clBufferPx
        memUsed[0] += nFrameOnGPU; // clBufferShiftX
        memUsed[0] += nFrameOnGPU; // clBufferShiftY
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGx
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGy
        memUsed[0] += width * height * magnification * magnification; // clBufferSRRF_AVG
        memUsed[0] += width * height * magnification * magnification; // clBufferSRRF_VAR

        // Memory on CPU ----
        memUsed[1] = 0;
        memUsed[1] += width * height * nSlices; // clBufferPx
        memUsed[1] += width * height * magnification * magnification; // clBufferSRRF // TODO: estimating RAM requirements properly

        memUsed[0] *= Float.SIZE / 8000000d;
        memUsed[1] *= Float.SIZE / 8000000d;

        return memUsed;
    }


}