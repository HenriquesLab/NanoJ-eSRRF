package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PointRoi;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import org.python.modules.math;

import java.awt.*;

import static java.lang.Math.min;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

public class liveSRRF_SARF_ implements PlugIn {

    // Basic formats
    private int magnification, nFrameSRRF, sensitivity, frameGap, nFrameOnGPU;
    private float fwhm, maxMemoryGPU;
    private boolean correctVibration, calculateAVG, calculateSTD, doFusion, getInterpolatedImage;
    private final int radiusCCM = 5;

    // Image formats
    private ImagePlus imp;
    private ImagePlus impCCM = null;

    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();


    @Override
    public void run(String arg) {

        // Get raw data
        imp = WindowManager.getCurrentImage(); // TODO: depending on the size of data and RAM, consider using Virtual Stack load
        if (imp == null) imp = IJ.openImage();
        imp.show();

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("LiveSRRF - Unfinished...");
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("FWHM (pixels, default: 2)", prefs.get("fwhm", 2), 2);
        gd.addNumericField("Sensitivity (default: 3)", prefs.get("fwhm", 3), 0);

        gd.addNumericField("# frames for SRRF (0 = auto)", prefs.get("nFrameSRRF", 0), 0);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("showAVGReconstruction", true));
        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("showSTDReconstruction", false));
        gd.addCheckbox("Fusion reconstruction (default: off)", prefs.get("showFusionReconstruction", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("showWFInterpolation", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addNumericField("Gap between SR frame (frames, default: 50)", prefs.get("frameGap_RA", 50), 0);

        gd.addMessage("-=-= Memory =-=-\n", headerFont);
        gd.addNumericField("Maximum amount of memory on GPU (MB, default: 1000)", prefs.get("maxMemGPU", 500), 2);

        MyDialogListener dl = new MyDialogListener(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        IJ.log("-------------------------------------");
        IJ.log("liveSRRF - SARF");
        IJ.log("Magnification: " + magnification);
        IJ.log("# frames for SRRF: " + nFrameSRRF);
        IJ.log("# frames gap: " + frameGap);
        IJ.log("# frames on GPU: " + nFrameOnGPU);
        IJ.log("# frame gap: " + predictMemoryUsed(nFrameOnGPU)[0] + " MB");


        // Get ready

        float[] shiftX = new float[nFrameOnGPU];
        float[] shiftY = new float[nFrameOnGPU];

        imp.setSlice(1);
        ImageProcessor ipRef = imp.getProcessor();
        ImageProcessor ipData;
//        ImagePlus impCCM = null;


        for (int s = 1; s <= nFrameOnGPU; s++) {
            // Check if user is cancelling calculation
            IJ.showProgress(s, nFrameOnGPU);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }

            // Grab the new frame from the list
            imp.setSlice(s);
            ipData = imp.getProcessor();

            // Estimate vibrations
            if (correctVibration) {
                int id = prof.startTimer();
                float[] shift = calculateShift(ipRef, ipData, radiusCCM, impCCM);
                shiftX[s] = shift[0];
                shiftY[s] = shift[1];

                System.out.println("Frame=" + s + " shiftX=" + shiftX[s] + " shiftY=" + shiftY[s]);
                prof.recordTime("Drift Estimation", prof.endTimer(id));
            }

        }


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

        calculateAVG = gd.getNextBoolean();
        calculateSTD = gd.getNextBoolean();
        doFusion = gd.getNextBoolean();
        getInterpolatedImage = gd.getNextBoolean();


        frameGap = (int) gd.getNextNumber();
        maxMemoryGPU = (int) gd.getNextNumber();

        //nBufferFrames = (int) ((float) nFrameSRRF / (float) frameGap);  // not used for SARF

        if (nFrameSRRF == 0) nFrameSRRF = imp.getImageStack().getSize();
        nFrameSRRF = min(imp.getImageStack().getSize(), nFrameSRRF);

        int maxnFrameOnGPU = 0;
        double[] memUsed = new double[2];
        while (memUsed[0] < (double) maxMemoryGPU) {
            maxnFrameOnGPU = maxnFrameOnGPU + 1;
            memUsed = predictMemoryUsed(maxnFrameOnGPU);
        }

        maxnFrameOnGPU = maxnFrameOnGPU - 1;
        if (maxnFrameOnGPU > 0) {
            goodToGo = true;
            nFrameOnGPU = (int) math.ceil((float) imp.getImageStack().getSize() / (float) maxnFrameOnGPU);
            nFrameOnGPU = (int) math.ceil(((float) imp.getImageStack().getSize() / (float) nFrameOnGPU));
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
        memUsed[0] += width * height * magnification * magnification; // clBufferRGC
        if (getInterpolatedImage) memUsed[0] += width * height * magnification * magnification; // clBufferInt

        // Memory on CPU ----
        memUsed[1] = 0;
        memUsed[1] += width * height * nSlices; // clBufferPx
        memUsed[1] += width * height * magnification * magnification; // clBufferSRRF // TODO: estimating RAM requirements properly

        memUsed[0] *= Float.SIZE / 8000000d;
        memUsed[1] *= Float.SIZE / 8000000d;

        return memUsed;
    }


    // --- Calculate shift using Cross-correlation matrix ---
    private float[] calculateShift(ImageProcessor ipRef, ImageProcessor ipData, int radius, ImagePlus impCCM) {

        FloatProcessor fpCCM = (FloatProcessor) calculateCrossCorrelationMap(ipRef, ipData, false);

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


}
