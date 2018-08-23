package nanoj.liveSRRF.gui;

import com.jogamp.opencl.CLDevice;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.PointRoi;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.liveSRRF.ErrorMapLiveSRRF;
import nanoj.liveSRRF.XYShiftCalculator;
import nanoj.liveSRRF.liveSRRF_CL;

import java.awt.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import static java.lang.Math.min;
import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

public class ParametersSweep_ implements PlugIn {

    // Basic formats
    private int magnification,
            nSlices,
            width,
            height,
            blockSize;

    private boolean calculateAVG,
            calculateSTD,
            showRecons,
            showErrorMaps,
            calculateFRC,
            calculateRSE,
            correctVibration;

    private float[] fwhmArray;
    private int[] sensitivityArray,
            nframeArray;

    private final String LiveSRRFVersion = "v0.4";
    private float[] shiftX, shiftY;

    // Image formats
    private ImagePlus imp;

    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();
    private liveSRRF_CL liveSRRF;

    public void run(String arg) {


        // Get raw data
        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        nSlices = imp.getImageStack().getSize(); // TODO: use nSlices to set maximum on #frames
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("liveSRRF - Parameters sweep " + LiveSRRFVersion);
        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));
        liveSRRF = new liveSRRF_CL();

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("liveSRRF - Parameters sweep " + LiveSRRFVersion);
        gd.addMessage("-=-= Fixed SRRF parameters =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addCheckbox("Correct vibration (default: off)", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Sweeping SRRF parameters =-=-\n", headerFont);
        gd.addMessage("Radius\n", headerFont);
        gd.addNumericField("Start", prefs.get("fwhm0", 2), 2);
        gd.addNumericField("Delta", prefs.get("deltafwhm", 0.5f), 2);
        gd.addNumericField("Number", prefs.get("n_fwhm", 5), 0);

        gd.addMessage("Sensitivity\n", headerFont);
        gd.addNumericField("Start", prefs.get("S0", 1), 0);
        gd.addNumericField("Delta", prefs.get("deltaS", 1), 0);
        gd.addNumericField("Number", prefs.get("n_S", 4), 0);

        gd.addMessage("# frames for SRRF\n", headerFont);
        gd.addNumericField("Start", prefs.get("nf0", 50), 0);
        gd.addNumericField("Delta", prefs.get("deltanf", 25), 0);
        gd.addNumericField("Number", prefs.get("n_nf", 3), 0);

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("calculateSTD", false));

        gd.addMessage("-=-= Output =-=-\n", headerFont);
        gd.addCheckbox("Show all reconstructions (default: on)", prefs.get("showRecons", true));
        gd.addCheckbox("Show all error maps (default: on)", prefs.get("showErrorMaps", true));
        gd.addCheckbox("Calculate FRC (default: off)", prefs.get("calculateFRC", false));
        gd.addCheckbox("Calculate RSE (default: off)", prefs.get("calculateRSE", false));

        gd.addMessage("Calculating FRC will split all dataset in two halves and therefore\n" +
                "the maximum number of frames will be half of the total frames\n" +
                "available in the dataset.");

        gd.addMessage("-=-= GPU processing =-=-\n", headerFont);
        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use\n" +
                "more resources and may slow down your computer.");

        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            liveSRRF.release();
            return;
        }

        grabSettings(gd);

        ImageStack imsAllRawData = imp.getImageStack();
        ImageStack imsThisRawData;

        ImageStack imsBuffer;
        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);
        ImageStack imsErrorMap = new ImageStack(width * magnification, height * magnification);

        ImageStack imsRMSE = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());


        int n_calculation = nframeArray.length * sensitivityArray.length * fwhmArray.length;

        IJ.log("Magnification: " + magnification);
        if (correctVibration) IJ.log("Vibration correction: on");
        else IJ.log("Vibration correction: off");

        IJ.log("Number of calculations planned: " + n_calculation);
        int r = 0;

        boolean userPressedEscape;
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF();
        FloatProcessor fpErrorMap;
        float[] pixelsRMSE;

        int maxnFrame = nframeArray[nframeArray.length - 1];
        shiftX = new float[maxnFrame];
        shiftY = new float[maxnFrame];
        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp, prof);

//        if (correctVibration) calculateShiftArray(1, maxnFrame);
        if (correctVibration) {
            shiftCalculator.calculateShiftArray(1, maxnFrame);
            shiftX = shiftCalculator.shiftX;
            shiftY = shiftCalculator.shiftY;
        }

        float[] shiftXtemp;
        float[] shiftYtemp;


        for (int nfi = 0; nfi < nframeArray.length; nfi++) {
            pixelsRMSE = new float[(fwhmArray.length) * (sensitivityArray.length)];

            shiftXtemp = new float[nframeArray[nfi]];
            shiftYtemp = new float[nframeArray[nfi]];

            for (int i = 0; i < nframeArray[nfi]; i++) {
                shiftXtemp[i] = shiftX[i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
                shiftYtemp[i] = shiftY[i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
            }

            for (int si = 0; si < sensitivityArray.length; si++) {
                for (int fi = 0; fi < fwhmArray.length; fi++) {

                    IJ.log("--------");
                    IJ.log("SRRF frame: " + r+1 + "/" + n_calculation);
                    IJ.showProgress(r, n_calculation);

                    // Check if user is cancelling calculation
                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        liveSRRF.release();
                        return;
                    }

                    IJ.log("Number of frame for SRRF: " + nframeArray[nfi]);
                    IJ.log("Radius: " + fwhmArray[fi] + " pixels");
                    IJ.log("Sensitivity: " + sensitivityArray[si]);

                    liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null);
                    liveSRRF.resetFramePosition();
                    liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

                    IJ.showStatus("Calculating SRRF image...");
                    for (int f = 1; f <= nframeArray[nfi]; f++) {
//                        imp.setSlice(f);
                        imsThisRawData = new ImageStack(width, height);
//                        imsThisRawData.addSlice(imp.getProcessor());
                        imsThisRawData.addSlice(imsAllRawData.getProcessor(f));
                        userPressedEscape = liveSRRF.calculateSRRF(imsThisRawData);

                        // Check if user is cancelling calculation
                        if (userPressedEscape) {
                            liveSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }
                    }

                    IJ.showStatus("Reading SRRF image...");
                    imsBuffer = liveSRRF.readSRRFbuffer();

                    IJ.showStatus("Optimising Sigma...");
                    errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(1), magnification);

                    IJ.showStatus("Calculating error map...");
                    fpErrorMap = errorMapCalculator.calculateErrorMap();
                    pixelsRMSE[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;

                    String label = "R=" + fwhmArray[fi] + "/S=" + sensitivityArray[si] + "/#fr=" + nframeArray[nfi];
                    if (showErrorMaps)
                        imsErrorMap.addSlice(label, fpErrorMap);
                    if (calculateAVG)
                        imsSRRFavg.addSlice(label, imsBuffer.getProcessor(1));
                    if (calculateSTD)
                        imsSRRFstd.addSlice(label, imsBuffer.getProcessor(2));
                    imsInt.addSlice(label, imsBuffer.getProcessor(3));

                    r++;
                }
            }
            imsRMSE.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSE), nfi + 1);
            imsRMSE.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);

        }

        // Release the GPU
        liveSRRF.release();

        //Display results
        if (calculateAVG) {
            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG)", imsSRRFavg);
            impSRRFavg.setCalibration(cal);
            IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
            impSRRFavg.show();
        }

        if (calculateSTD) {
            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD)", imsSRRFstd);
            impSRRFstd.setCalibration(cal);
            IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
            impSRRFstd.show();
        }


        ImagePlus impInt = new ImagePlus(imp.getTitle() + " - Interpolated image", imsInt);
        impInt.setCalibration(cal);
        IJ.run(impInt, "Enhance Contrast", "saturated=0.5");
        impInt.show();

        if (showErrorMaps) {
            ImagePlus impErrorMap = new ImagePlus(imp.getTitle() + " - ErrorMap", imsErrorMap);
            impErrorMap.setCalibration(cal);
            IJ.run(impErrorMap, "Enhance Contrast", "saturated=0.5");
            applyLUT_SQUIRREL_Errors(impErrorMap);
            impErrorMap.show();
        }

        if (calculateRSE) {
            ImagePlus impRMSE = new ImagePlus(imp.getTitle() + " - RMSE sweep map", imsRMSE);
            impRMSE.setCalibration(cal);
            IJ.run(impRMSE, "Enhance Contrast", "saturated=0.5");
            applyLUT_SQUIRREL_Errors(impRMSE);
            impRMSE.show();
        }

        IJ.log("-------------------------------------");
        IJ.log("RAM used: " + IJ.freeMemory());
        IJ.log("Bye-bye !");

        // Run garbage collector
        System.gc();

    }

    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    //    --- Grab settings ---
    private void grabSettings(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();
        correctVibration = gd.getNextBoolean();

        float fwhm0 = (float) gd.getNextNumber();
        float deltafwhm = (float) gd.getNextNumber();
        int n_fwhm = (int) gd.getNextNumber();

        int S0 = (int) gd.getNextNumber();
        int deltaS = (int) gd.getNextNumber();
        int n_S = (int) gd.getNextNumber();

        int nf0 = (int) gd.getNextNumber();
        int deltanf = (int) gd.getNextNumber();
        int n_nf = (int) gd.getNextNumber();

        calculateAVG = gd.getNextBoolean();
        calculateSTD = gd.getNextBoolean();

        showRecons = gd.getNextBoolean();
        showErrorMaps = gd.getNextBoolean();

        calculateFRC = gd.getNextBoolean();
        calculateRSE = gd.getNextBoolean();

        blockSize = (int) gd.getNextNumber();

        fwhmArray = new float[n_fwhm];
        for (int i = 0; i < n_fwhm; i++) {
            fwhmArray[i] = fwhm0 + i * deltafwhm;
        }

        sensitivityArray = new int[n_S];
        for (int i = 0; i < n_S; i++) {
            sensitivityArray[i] = S0 + i * deltaS;
        }

        nframeArray = new int[n_nf];
        for (int i = 0; i < n_nf; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }

        prefs.set("magnification", magnification);
        prefs.set("correctVibration", correctVibration);

        prefs.set("fwhm0", fwhm0);
        prefs.set("deltafwhm", deltafwhm);
        prefs.set("n_fwhm", n_fwhm);

        prefs.set("S0", S0);
        prefs.set("deltaS", deltaS);
        prefs.set("n_S", n_S);

        prefs.set("nf0", nf0);
        prefs.set("deltanf", deltanf);
        prefs.set("n_nf", n_nf);

        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateSTD", calculateSTD);

        prefs.set("showRecons", showRecons);
        prefs.set("showErrorMaps", showErrorMaps);
        prefs.set("calculateFRC", calculateFRC);
        prefs.set("calculateRSE", calculateRSE);


        prefs.set("blockSize", blockSize);


        prefs.save();
    }


}
