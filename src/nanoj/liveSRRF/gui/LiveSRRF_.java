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
import nanoj.liveSRRF.SRRF2CL;

import java.awt.*;

import static java.lang.Math.*;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;
import static nanoj.liveSRRF.SRRF2CL.BIN_4;
import static nanoj.liveSRRF.SRRF2CL.predictMemoryUsed;

public class LiveSRRF_ implements PlugIn {

    private Font headerFont = new Font("Arial", Font.BOLD, 16);

    private String user = "HenriquesLab";
    //private String user = "AudreySalles";
    //private String user = "PaulReynolds";
    //private String user = "GrahamDellaire";
    private String version = "20180418-"+user;

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJUsageTracker tracker = new NanoJUsageTracker("NanoJ-LiveSRRF", version, "UA-61590656-4");
    private NanoJProfiler prof = new NanoJProfiler();

    private ImageStack imsSRRF;
    private ImagePlus impCCM = null;

    private ImagePlus imp;
    private int magnification, nFrames, nTimeLags;
    private float fwhm;
    private boolean correctVibration, correctSCMOS, showAllReconstructions, doBin2, doBin4;

    @Override
    public void run(String arg) {
        tracker.logUsage("LiveSRRF_");

        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Live SRRF (aka SRRF2) - In Development (for "+user+")");
        gd.addNumericField("Magnification (default: 4)", prefs.get("magnification", 4), 0);
        gd.addNumericField("FWHM (pixels, default: 3)", prefs.get("fwhm", 3), 2);
        gd.addNumericField("Frames_per_time-point (0 - auto)", prefs.get("nFrames", 0), 0);
        gd.addCheckbox("Correct vibration", prefs.get("correctVibration", false));
        gd.addCheckbox("Correct sCMOS patterning", prefs.get("correctSCMOS", false));
        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("Show_all reconstruction", prefs.get("showAllReconstructions", false));
        gd.addNumericField("Accumulative_timelags (default: 5)", prefs.get("nTimeLags", 5), 0);
        gd.addCheckbox("Calculate for bin 2 (default: off)", prefs.get("doBin2", false));
        gd.addCheckbox("Calculate for bin 4 (default: off)", prefs.get("doBin4", false));
        gd.addMessage("-=-= Advice =-=-\n", headerFont);
        gd.addMessage(
                "SRRF2 is a GPU resources-hog, if it fails to run consider\n" +
                     "reducing the magnification. For larger temporal datasets,\n" +
                     "consider doing batches of 100 frames at a time, then do a \n" +
                     "a mean projecting for them.");

        gd.addMessage("-=-= CONFIDENTIAL =-=-\n", headerFont);
        gd.addMessage(
                "This is a preview of the 'in development' version of the \n" +
                        "SRRF2 engine developed by the Henriques lab @ UCL.\n" +
                        "It is only meant to be used by researchers who received\n" +
                        "a direct email by Ricardo Henriques.");

        MyDialogListener dl = new MyDialogListener(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        gd.showDialog();
        if (gd.wasCanceled()) return;
        //grabSettings(gd);

        if (correctSCMOS) {
            IJ.showMessage("SCMOS correction is disabled for now on purpose =)");
        }

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();
        int wM = ims.getWidth() * magnification;
        int hM = ims.getHeight() * magnification;
        int nPixelsM = wM * hM;

        ImagePlus impSRRF = new ImagePlus(imp.getTitle()+" - SRRF2");
        impSRRF.copyScale(imp); // make sure we copy the pixel sizes correctly accross
        Calibration cal = impSRRF.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;

        imsSRRF = new ImageStack(wM, hM);
        //impSRRF.setStack(imsSRRF);
        ImageStack imsRawDataBuffer = new ImageStack(w, h);
        boolean firstTime = true;

        ImageProcessor ipRef = null; // reference slide for Cross-Correlation and vibration correction

        SRRF2CL srrf2CL = new SRRF2CL(w, h, nFrames, magnification, fwhm, nTimeLags, BIN_4);

        float[] shiftX = new float[nFrames];
        float[] shiftY = new float[nFrames];
        int counter = 0;
        int nReconstructions = 1;

        //////////////////////////////////////
        // !!! MAIN LOOP THROUGH FRAMES !!! //
        //////////////////////////////////////

        for (int s=1; s<=nSlices; s++) {
            // Check if user is cancelling calculation
            IJ.showProgress(s, nSlices);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                srrf2CL.release();
                return;
            }

            // Grab the new frame from the list
            imp.setSlice(s);
            ImageProcessor ip = imp.getProcessor();
            imsRawDataBuffer.addSlice(ip);

            // Estimate vibrations
            if (correctVibration) {
                System.out.println("New reference..."+counter);
                int id = prof.startTimer();
                if (counter == 0) {
                    ipRef = ip.duplicate();
                    shiftX[counter] = 0;
                    shiftY[counter] = 0;
                }
                else {
                    float[] shift = calculateShift(ipRef, ip, 5);
                    shiftX[counter] = shift[0];
                    shiftY[counter] = shift[1];
                }
                System.out.println("Frame="+s+" shiftX="+shiftX[counter]+" shiftY="+shiftY[counter]);
                prof.recordTime("Drift Estimation", prof.endTimer(id));
            }

            if (counter == nFrames-1 || s == nSlices) {
                int id = prof.startTimer();
                ImageStack imsResults = srrf2CL.calculateSRRF(imsRawDataBuffer, shiftX, shiftY);
                nReconstructions = imsResults.getSize();
                if (!showAllReconstructions) {
                    imsSRRF.addSlice(imsResults.getProcessor(nReconstructions));
                    imsSRRF.setSliceLabel(srrf2CL.reconstructionLabel.get(nReconstructions-1), imsSRRF.getSize());
                }
                else {
                    for (int r = 1; r <= imsResults.getSize(); r++) {
                        imsSRRF.addSlice(imsResults.getProcessor(r));
                        imsSRRF.setSliceLabel(srrf2CL.reconstructionLabel.get(r - 1), imsSRRF.getSize());
                    }
                }
                if (firstTime) {
                    impSRRF.setStack(imsSRRF);
                    impSRRF.show();
                    impSRRF.setSlice(imsSRRF.getSize());
                    IJ.run(impSRRF, "Enhance Contrast", "saturated=0.35");
                    firstTime = false;
                }
                else {
                    impSRRF.setSlice(imsSRRF.getSize());
                }

                // reset buffers
                imsRawDataBuffer = new ImageStack(w, h);
                counter = 0;
                prof.recordTime("full SRRF-frame calculation", prof.endTimer(id));
            }
            else counter++;
        }

        srrf2CL.release(); // Release the GPU!!!
        IJ.log(prof.report());

        // Show final rendering...
        impSRRF.setStack(imsSRRF);
        IJ.run(impSRRF, "Enhance Contrast", "saturated=0.5");
        impSRRF.setTitle(imp.getTitle()+" - SRRF2");

        if (showAllReconstructions) {
            int nSRRFFrames = imsSRRF.getSize() / nReconstructions;
            IJ.run(impSRRF, "Stack to Hyperstack...", "order=xyczt(default) channels="+nReconstructions+" slices=1 frames="+nSRRFFrames+" display=Grayscale");
        }
    }

    private boolean grabSettings(GenericDialog gd) {
        magnification = (int) gd.getNextNumber();
        fwhm = (float) gd.getNextNumber();
        nFrames = (int) gd.getNextNumber();
        correctVibration = gd.getNextBoolean();
        correctSCMOS = gd.getNextBoolean();
        showAllReconstructions = gd.getNextBoolean();
        nTimeLags = (int) gd.getNextNumber();
        doBin2 = gd.getNextBoolean();
        doBin4 = gd.getNextBoolean();

        if (nTimeLags<1) return false;

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrames", nFrames);
        prefs.set("correctVibration", correctVibration);
        prefs.set("correctSCMOS", correctSCMOS);

        prefs.set("showAllReconstructions", showAllReconstructions);
        prefs.set("nTimeLags", nTimeLags);
        prefs.set("doBin2", doBin2);
        prefs.set("doBin4", doBin4);
        prefs.save();

        if (nFrames == 0) nFrames = imp.getImageStack().getSize();
        nFrames = min(imp.getImageStack().getSize(), nFrames);

        return true;
    }

    class MyDialogListener implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {

            boolean goodToGo = grabSettings(gd);
            ImageStack ims = imp.getImageStack();

            double memUsed = predictMemoryUsed(ims.getWidth(), ims.getHeight(), nFrames, magnification, nTimeLags, BIN_4);
            IJ.showStatus("SRRF2 - Predicted used GPU memory: "+Math.round(memUsed)+"MB");

            return goodToGo;
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


}