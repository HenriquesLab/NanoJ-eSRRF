package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
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

public class LiveSRRF2_ implements PlugIn {

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
    private boolean showAllReconstructions;

    @Override
    public void run(String arg) {
        tracker.logUsage("LiveSRRF_");

        ImagePlus imp = WindowManager.getCurrentImage();
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

        gd.addMessage("-=-= CONFIDENTIAL =-=-\n", headerFont);
        gd.addMessage(
                "This is a preview of the 'in development' version of the \n" +
                        "SRRF2 engine developed by the Henriques lab @ UCL.\n" +
                        "It is only meant to be used by researchers who received\n" +
                        "a direct email by Ricardo Henriques.");

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int magnification = (int) gd.getNextNumber();
        float fwhm = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        boolean correctVibration = gd.getNextBoolean();
        boolean correctSCMOS = gd.getNextBoolean();
        boolean showAllReconstructions = gd.getNextBoolean();

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrames", nFrames);
        prefs.set("correctVibration", correctVibration);
        prefs.set("correctSCMOS", correctSCMOS);

        prefs.set("showAllReconstructions", showAllReconstructions);
        prefs.save();

        if (correctSCMOS) {
            IJ.showMessage("SCMOS correction is disabled for now on purpose =)");
        }

        if (nFrames == 0) nFrames = imp.getImageStack().getSize();
        nFrames = min(imp.getImageStack().getSize(), nFrames);

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

        SRRF2CL srrf2CL = new SRRF2CL(w, h, nFrames, magnification, fwhm);

        float[] shiftX = new float[nFrames];
        float[] shiftY = new float[nFrames];
        int counter = 0;

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
                int nReconstructions = imsResults.getSize();
                if (!showAllReconstructions) {
                    imsSRRF.addSlice(imsResults.getProcessor(nReconstructions));
                    imsSRRF.setSliceLabel(srrf2CL.reconstructionLabel[nReconstructions-1], imsSRRF.getSize());
                }
                else {
                    for (int r = 1; r <= imsResults.getSize(); r++) {
                        imsSRRF.addSlice(imsResults.getProcessor(r));
                        imsSRRF.setSliceLabel(srrf2CL.reconstructionLabel[r - 1], imsSRRF.getSize());
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

        impSRRF.setStack(imsSRRF);
        IJ.run(impSRRF, "Enhance Contrast", "saturated=0.35");
        impSRRF.setTitle(imp.getTitle()+" - SRRF2");

        srrf2CL.release();
        IJ.log(prof.report());
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