package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core.java.image.analysis.FRC;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.liveSRRF.ErrorMapLiveSRRF;
import nanoj.liveSRRF.XYShiftCalculator;
import nanoj.liveSRRF.liveSRRF_CL;

import java.awt.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;

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
            showRSC,
            calculateFRC,
            calculateRSE,
            calculateRSP,
            fixSigma,
            correctVibration;

    private float[] fwhmArray;
    private int[] sensitivityArray,
            nframeArray;

    private float fixedSigma;

    private final String LiveSRRFVersion = "v0.7";
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

        nSlices = imp.getImageStack().getSize();
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
        gd.addToSameRow();
        gd.addNumericField("Delta", prefs.get("deltafwhm", 0.5f), 2);
        gd.addToSameRow();
        gd.addNumericField("Number", prefs.get("n_fwhm", 5), 0);

        gd.addMessage("Sensitivity\n", headerFont);
        gd.addNumericField("Start", prefs.get("S0", 1), 0);
        gd.addToSameRow();
        gd.addNumericField("Delta", prefs.get("deltaS", 1), 0);
        gd.addToSameRow();
        gd.addNumericField("Number", prefs.get("n_S", 4), 0);

        gd.addMessage("# frames for SRRF\n", headerFont);
        gd.addNumericField("Start", prefs.get("nf0", 50), 0);
        gd.addToSameRow();
        gd.addNumericField("Delta", prefs.get("deltanf", 25), 0);
        gd.addToSameRow();
        gd.addNumericField("Number", prefs.get("n_nf", 3), 0);

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addToSameRow();
        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("calculateSTD", false));

        gd.addMessage("-=-= Output =-=-\n", headerFont);
        gd.addCheckbox("Show all reconstructions (default: on)", prefs.get("showRecons", true));
        gd.addCheckbox("Show all error maps (default: on)", prefs.get("showErrorMaps", true));
        gd.addCheckbox("Show all rescaled reconstructions (default: off)", prefs.get("showRSC", false));

        gd.addCheckbox("Calculate RSE (default: off)", prefs.get("calculateRSE", false));
        gd.addToSameRow();
        gd.addCheckbox("Calculate RSP (default: off)", prefs.get("calculateRSP", false));

        gd.addCheckbox("Fix sigma (default: off)", prefs.get("fixSigma", false));
        gd.addToSameRow();
        gd.addNumericField("Sigma", prefs.get("fixedSigma", 5), 2);

        gd.addCheckbox("Calculate FRC (default: off)", prefs.get("calculateFRC", false));

        gd.addMessage("Calculating FRC will split all dataset in two halves and therefore\n" +
                "the maximum number of frames will be half of the total frames\n" +
                "available in the dataset.");

        gd.addMessage("-=-= GPU processing =-=-\n", headerFont);
        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use\n" +
                "more resources and may slow down your computer.");

        gd.addHelp("https://www.youtube.com/watch?v=Vs6awg-BJHo"); // If you don't know how to sweep

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
        ImageStack imsRSC = new ImageStack(width * magnification, height * magnification);


        ImageStack imsRMSE = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCC = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);

        ImageStack imsFRCresolutionAVG = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsFRCresolutionSTD = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);


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

        if (fixSigma) IJ.log("Sigma is fixed to " + fixedSigma + " pixels.");
        else IJ.log("Sigma is optimised for each reconstructions.");

        IJ.log("Number of calculations planned: " + n_calculation);
        int r = 0;

        boolean userPressedEscape;
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF(imp, magnification, fixSigma, fixedSigma);
        FloatProcessor fpErrorMap;
        float[] pixelsRMSE;
        float[] pixelsPPMCC;

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

        FRC frcCalculator = new FRC();
        FloatProcessor ipOddAVG;
        FloatProcessor ipEvenAVG;
        FloatProcessor ipOddSTD;
        FloatProcessor ipEvenSTD;

        float[] pixelsFRCresolutionAVG;
        float[] pixelsFRCresolutionSTD;



        for (int nfi = 0; nfi < nframeArray.length; nfi++) {
            pixelsRMSE = new float[(fwhmArray.length) * (sensitivityArray.length)];
            pixelsPPMCC = new float[(fwhmArray.length) * (sensitivityArray.length)];

            shiftXtemp = new float[nframeArray[nfi]];
            shiftYtemp = new float[nframeArray[nfi]];

            pixelsFRCresolutionAVG = new float[(fwhmArray.length) * (sensitivityArray.length)];
            pixelsFRCresolutionSTD = new float[(fwhmArray.length) * (sensitivityArray.length)];

            for (int i = 0; i < nframeArray[nfi]; i++) {
                shiftXtemp[i] = shiftX[i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
                shiftYtemp[i] = shiftY[i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
            }

            for (int si = 0; si < sensitivityArray.length; si++) {
                for (int fi = 0; fi < fwhmArray.length; fi++) {

                    IJ.log("--------");
                    IJ.log("SRRF frame: " + (r + 1) + "/" + n_calculation);
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

//                    liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true);
//                    liveSRRF.resetFramePosition();
//                    liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

//                    IJ.showStatus("Calculating SRRF image...");
//                    for (int f = 1; f <= nframeArray[nfi]; f++) {
////                        imp.setSlice(f);
//                        imsThisRawData = new ImageStack(width, height);
////                        imsThisRawData.addSlice(imp.getProcessor());
//                        imsThisRawData.addSlice(imsAllRawData.getProcessor(f));
//                        userPressedEscape = liveSRRF.calculateSRRF(imsThisRawData);
//
//                        // Check if user is cancelling calculation
//                        if (userPressedEscape) {
//                            liveSRRF.release();
//                            IJ.log("-------------------------------------");
//                            IJ.log("Reconstruction aborted by user.");
//                            return;
//                        }
//                    }
//
//                    IJ.showStatus("Reading SRRF image...");
////                    imsBuffer = liveSRRF.readSRRFbuffer();
//                    liveSRRF.readSRRFbuffer();

                    if (calculateFRC) {

                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 1);
                        imsBuffer = liveSRRF.imsSRRF;

                        ipOddAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                        ipOddSTD = imsBuffer.getProcessor(2).convertToFloatProcessor();

                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 2);
                        imsBuffer = liveSRRF.imsSRRF;

                        ipEvenAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                        ipEvenSTD = imsBuffer.getProcessor(2).convertToFloatProcessor();

                        pixelsFRCresolutionAVG[fwhmArray.length * si + fi] = (float) frcCalculator.calculateFireNumber(ipOddAVG, ipEvenAVG, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        pixelsFRCresolutionSTD[fwhmArray.length * si + fi] = (float) frcCalculator.calculateFireNumber(ipOddSTD, ipEvenSTD, FRC.ThresholdMethod.FIXED_1_OVER_7);

                    } else {

                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 0);
                        imsBuffer = liveSRRF.imsSRRF;

                    }


                    IJ.showStatus("Optimising Sigma...");
                    errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(1)); // TODO: currently only doing it on AVG images

                    IJ.showStatus("Calculating error map...");
                    fpErrorMap = errorMapCalculator.calculateErrorMap();
                    pixelsRMSE[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                    pixelsPPMCC[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                    String label = "R=" + fwhmArray[fi] + "/S=" + sensitivityArray[si] + "/#fr=" + nframeArray[nfi];
                    if (showErrorMaps) imsErrorMap.addSlice(label, fpErrorMap);
                    if (showErrorMaps) imsRSC.addSlice(label, errorMapCalculator.fpSRC);

                    if (calculateAVG) imsSRRFavg.addSlice(label, imsBuffer.getProcessor(1));
                    if (calculateSTD) imsSRRFstd.addSlice(label, imsBuffer.getProcessor(2));
                    imsInt.addSlice(label, imsBuffer.getProcessor(3));

                    r++;
                }
            }

            if (calculateRSE) {
                imsRMSE.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSE), nfi + 1);
                imsRMSE.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
            }

            if (calculateRSP) {
                imsPPMCC.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCC), nfi + 1);
                imsPPMCC.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
            }

            if (calculateFRC) {
                imsFRCresolutionAVG.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionAVG), nfi + 1);
                imsFRCresolutionAVG.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
            }

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

        if (showRSC) {
            ImagePlus impRSC = new ImagePlus(imp.getTitle() + " - rescaled liveSRRF", imsRSC);
            impRSC.setCalibration(cal);
            IJ.run(impRSC, "Enhance Contrast", "saturated=0.5");
            impRSC.show();
        }

        if (calculateRSE) {
            ImagePlus impRMSE = new ImagePlus(imp.getTitle() + " - RMSE sweep map", imsRMSE);
            impRMSE.setCalibration(cal);
            IJ.run(impRMSE, "Enhance Contrast", "saturated=0.5");
            applyLUT_SQUIRREL_Errors(impRMSE);
            impRMSE.show();
        }

        if (calculateRSP) {
            ImagePlus impPPMCC = new ImagePlus(imp.getTitle() + " - PMC sweep map", imsPPMCC);
            impPPMCC.setCalibration(cal);
            IJ.run(impPPMCC, "Enhance Contrast", "saturated=0.5");
            applyLUT_SQUIRREL_Errors(impPPMCC);
            impPPMCC.show();
        }

        if (calculateFRC) {
            ImagePlus impFRCresolution = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map", imsFRCresolutionAVG);
            impFRCresolution.setCalibration(cal);
            IJ.run(impFRCresolution, "Enhance Contrast", "saturated=0.5");
            applyLUT_SQUIRREL_Errors(impFRCresolution);
            impFRCresolution.show();
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
        showRSC = gd.getNextBoolean();

        calculateRSE = gd.getNextBoolean();
        calculateRSP = gd.getNextBoolean();

        fixSigma = gd.getNextBoolean();
        fixedSigma = (float) gd.getNextNumber();

        calculateFRC = gd.getNextBoolean();

        blockSize = (int) gd.getNextNumber();

        fwhmArray = new float[n_fwhm];
        for (int i = 0; i < n_fwhm; i++) {
            fwhmArray[i] = fwhm0 + i * deltafwhm;
        }

        sensitivityArray = new int[n_S];
        for (int i = 0; i < n_S; i++) {
            sensitivityArray[i] = S0 + i * deltaS;
        }

        // Check that nf does not exceed nSlices
        int n_nfToUse;
        if (calculateFRC) n_nfToUse = Math.min((nSlices / 2 - nf0) / deltanf + 1, n_nf); // TODO: review the +1
        else n_nfToUse = Math.min((nSlices - nf0) / deltanf + 1, n_nf);
        // TODO: ABORT IF NOT enough frames in original stack for even 1 analysis

        nframeArray = new int[n_nfToUse];
        for (int i = 0; i < n_nfToUse; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }

        //        ArrayList<Integer> nframeListArray = new ArrayList<Integer>();
//        int i = 0;
//        int nf = nf0;
//        while (nf <= nSlices && i < n_nf){
//            nframeListArray.add(nf);
//            nf = nf0 + i * deltanf;
//            i++;
//        }

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
        prefs.set("n_nf", n_nf); // save the one set by the user, not the calculated one

        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateSTD", calculateSTD);

        prefs.set("showRecons", showRecons);
        prefs.set("showErrorMaps", showErrorMaps);
        prefs.set("showRSC", showRSC);

        prefs.set("calculateRSE", calculateRSE);
        prefs.set("calculateRSP", calculateRSP);

        prefs.set("fixSigma", fixSigma);
        prefs.set("fixedSigma", fixedSigma);

        prefs.set("calculateFRC", calculateFRC);

        prefs.set("blockSize", blockSize);


        prefs.save();
    }


    public void calculateLiveSRRFsingleframeLoad(ImageStack imsAllRawData, int nf, int mode) {

        ImageStack imsThisRawData;
        boolean userPressedEscape;
        int fmode;

        IJ.showStatus("Calculating SRRF image...");
        for (int f = 1; f <= nframeArray[nf]; f++) {
            imsThisRawData = new ImageStack(width, height);

            if (mode == 0) fmode = f;  // no FRC
            else if (mode == 1) fmode = 2 * (f - 1) + 1; // FRC odd frames
            else fmode = 2 * f; // FRC even frames

            imsThisRawData.addSlice(imsAllRawData.getProcessor(fmode));
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
        liveSRRF.readSRRFbuffer();

    }


}
