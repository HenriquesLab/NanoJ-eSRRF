package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.DialogListener;
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
import nanoj.liveSRRF.LiveSRRF_CL;

import java.awt.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_Errors;
import static nanoj.core.java.tools.NJ_LUT.applyLUT_SQUIRREL_FRC;


public class ParametersSweep_ implements PlugIn {

    // Basic formats
    private int magnification,
            nSlices,
            width,
            height,
            blockSize,
            cropSize;

    private boolean showRecons,
            calculateAVG,
            calculateVAR,
            calculateTAC2,
            getInterpolatedImage,
            showErrorMaps,
            showRSC,
            calculateFRC,
            calculateRSE,
            calculateRSP,
            doErrorMapping,
            fixSigma,
            correctVibration,
            cropBorder,
            singleFrameLoad,
            doFHTinterpolation,
            previousAdvSettings = false;

    private final boolean DEBUG = false;

    private float[] fwhmArray;
    private int[] sensitivityArray,
            nframeArray;

    private float fixedSigma;

    private final String eSRRFVersion = "v1.5.0";
    private float[] shiftX, shiftY;

    private String imageTitle;


    // Advanced formats
    private final NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private final NanoJProfiler prof = new NanoJProfiler();

    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        imageTitle = imp.getTitle();

        nSlices = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();


        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("eSRRF - Parameters sweep " + eSRRFVersion);
        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));

        // Initialise the eSRRF engine (creates a context really)
        LiveSRRF_CL eSRRF = new LiveSRRF_CL(DEBUG);

        // Set advanced variables to default or preset values
        showErrorMaps = prefs.get("showErrorMaps", true);
        showRSC = prefs.get("showRSC", false);
        calculateRSE = prefs.get("calculateRSE", false);
        calculateRSP = prefs.get("calculateRSP", false);
        fixSigma = prefs.get("fixSigma", false);
        fixedSigma = prefs.get("fixedSigma", 1);
        cropBorder = prefs.get("cropBorder", true);
        cropSize = (int) prefs.get("cropSize", 3);
        calculateFRC = prefs.get("calculateFRC", false);
        doFHTinterpolation = prefs.get("doFHTinterpolation", true);

        blockSize = (int) prefs.get("blockSize", 20000);
        singleFrameLoad = prefs.get("singleFrameLoad", true);


        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("eSRRF - Parameters sweep " + eSRRFVersion);
        gd.addMessage("-=-= eSRRF reconstruction =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("VAR reconstruction (default: off)", prefs.get("calculateVAR", false));
        gd.addCheckbox("TAC2 reconstruction (default: off)", prefs.get("calculateTAC2", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("getInterpolatedImage", false));

        gd.addCheckbox("Correct vibration (default: off)", prefs.get("correctVibration", false));
        gd.addCheckbox("Show all reconstructions (default: on)", prefs.get("showRecons", true));


        gd.addMessage("-=-= Sweeping eSRRF parameters =-=-\n", headerFont);
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

        gd.addMessage("-=-= Advanced settings =-=-\n", headerFont);
        gd.addCheckbox("Show advanced settings", false);

        MyDialogListenerMainGUI dl = new MyDialogListenerMainGUI(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);

        gd.addHelp("https://www.youtube.com/watch?v=Vs6awg-BJHo"); // If you don't know how to sweep
        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            eSRRF.release();
            return;
        }

        boolean goodToGo = grabSettings(gd);
        if (!goodToGo){
            eSRRF.release();
            IJ.log("-------------------------------------");
            IJ.log("Aborted: Not enough frames in the stack for this analysis.");
            return;
        }

        // Log the settings used
        IJ.log("Magnification: " + magnification);
        if (correctVibration) IJ.log("Vibration correction: on");
        else IJ.log("Vibration correction: off");

        if (doErrorMapping) {
            IJ.log("Error mapping paramaters:");
            if (fixSigma) IJ.log("Sigma is fixed to " + fixedSigma + " pixels.");
            else IJ.log("Sigma is optimised for each reconstructions.");
            if (cropBorder) IJ.log("The borders will be cropped by "+cropSize+" pixels");

        }

        int n_calculation = nframeArray.length * sensitivityArray.length * fwhmArray.length;
        IJ.log("Number of calculations planned: " + n_calculation);
        if (singleFrameLoad) IJ.log("Using single frame loading...");

        // Get the calibration from the raw data and set the calibration for the reconstructions (needs to happen after setting magnification)
        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

        // Initialising the Main output variables
        ImageStack imsAllRawData = imp.getImageStack();
        ImageStack imsBuffer;
        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFvar = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFtac2 = new ImageStack(width * magnification, height * magnification);

        // Initialising the Error maps variables
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF(imp, magnification, fixSigma, fixedSigma, cropBorder, cropSize);

        float[] pixelsRMSEavg = null;
        float[] pixelsPPMCCavg = null;
        float[] pixelsRMSEvar = null;
        float[] pixelsPPMCCvar = null;
        float[] pixelsRMSEtac2 = null;
        float[] pixelsPPMCCtac2 = null;

        int widthSquirrel;
        int heightSquirrel;

        if (cropBorder){
            widthSquirrel = (width - 2*cropSize) * magnification;
            heightSquirrel = (height - 2*cropSize) * magnification;
        }
        else{
            widthSquirrel = width * magnification;
            heightSquirrel = height * magnification;
        }


        ImageStack imsErrorMapAVG = new ImageStack(widthSquirrel, heightSquirrel);
        ImageStack imsErrorMapVAR = new ImageStack(widthSquirrel, heightSquirrel);
        ImageStack imsErrorMapTAC2 = new ImageStack(widthSquirrel, heightSquirrel);

        ImageStack imsRSCavg = new ImageStack(widthSquirrel, heightSquirrel);
        ImageStack imsRSCvar = new ImageStack(widthSquirrel, heightSquirrel);
        ImageStack imsRSCtac2 = new ImageStack(widthSquirrel, heightSquirrel);

        // These images are small so it doesn't matter if they are initialised in any case
        ImageStack imsRMSEavg = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCCavg = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsRMSEvar = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCCvar = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsRMSEtac2 = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCCtac2 = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);


        // Initialising the shift correction variables
        int maxnFrame;
        if (calculateFRC) maxnFrame = 2 * nframeArray[nframeArray.length - 1];
        else maxnFrame = nframeArray[nframeArray.length - 1];

        shiftX = new float[maxnFrame];
        shiftY = new float[maxnFrame];
        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp, prof);

        if (correctVibration) {
            shiftCalculator.calculateShiftArray(1, maxnFrame);
            shiftX = shiftCalculator.shiftX;
            shiftY = shiftCalculator.shiftY;
        }

        float[] shiftXtemp;
        float[] shiftYtemp;

        float[] shiftXtempOdd;
        float[] shiftYtempOdd;

        float[] shiftXtempEven;
        float[] shiftYtempEven;


        // Initialising the FRC variab;es
        FRC frcCalculator = new FRC();
        FloatProcessor fpOddAVG;
        FloatProcessor fpEvenAVG;
        FloatProcessor fpOddVAR;
        FloatProcessor fpEvenVAR;
        FloatProcessor fpOddTAC2;
        FloatProcessor fpEvenTAC2;

        float[] pixelsFRCresolutionAVG = null;
        float[] pixelsFRCresolutionVAR = null;
        float[] pixelsFRCresolutionTAC2 = null;

        // These images are small so it doesn't matter if they are initialised in any case
        ImageStack imsFRCresolutionAVG = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsFRCresolutionVAR = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsFRCresolutionTAC2 = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);

        // Set the number of reconstuctions
        int r = 0;
        boolean userPressedEscape;

        for (int nfi = 0; nfi < nframeArray.length; nfi++) {

            if (doErrorMapping) {
                pixelsRMSEavg = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsPPMCCavg = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsRMSEvar = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsPPMCCvar = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsRMSEtac2 = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsPPMCCtac2 = new float[(fwhmArray.length) * (sensitivityArray.length)];
            }

            if (calculateFRC) {
                pixelsFRCresolutionAVG = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsFRCresolutionVAR = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsFRCresolutionTAC2 = new float[(fwhmArray.length) * (sensitivityArray.length)];
            }


            shiftXtemp = new float[nframeArray[nfi]];
            shiftYtemp = new float[nframeArray[nfi]];

            shiftXtempOdd = new float[nframeArray[nfi]];
            shiftYtempOdd = new float[nframeArray[nfi]];
            shiftXtempEven = new float[nframeArray[nfi]];
            shiftYtempEven = new float[nframeArray[nfi]];

            // Re-adjust the array of shifts depending on the frames used
            if (calculateFRC) {

                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXtempOdd[i] = shiftX[2 * i + 1];
                    shiftYtempOdd[i] = shiftY[2 * i + 1];
                }

                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXtempEven[i] = shiftX[2 * i];
                    shiftYtempEven[i] = shiftY[2 * i];
                }

            } else {
                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXtemp[i] = shiftX[i];
                    shiftYtemp[i] = shiftY[i];
                }
            }

            int nFramesOnGPU;
            if (singleFrameLoad) nFramesOnGPU = 1;
            else nFramesOnGPU = nframeArray[nfi];


            for (int si = 0; si < sensitivityArray.length; si++) {
                for (int fi = 0; fi < fwhmArray.length; fi++) {

                    IJ.log("--------");
                    IJ.log("SRRF frame: " + (r + 1) + "/" + n_calculation);
                    IJ.showProgress(r, n_calculation);

                    // Check if user is cancelling calculation
                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        eSRRF.release();
                        IJ.log("-------------------------------------");
                        IJ.log("Reconstruction aborted by user.");
                        return;
                    }

                    IJ.log("Number of frame for SRRF: " + nframeArray[nfi]);
                    IJ.log("Radius: " + fwhmArray[fi] + " pixels");
                    IJ.log("Sensitivity: " + sensitivityArray[si]);

                    String label = "R=" + fwhmArray[fi] + "/S=" + sensitivityArray[si] + "/#fr=" + nframeArray[nfi];


                    // FRC resolution estimation
                    if (calculateFRC) {

                        // Calculate and get the reconstruction from the odd frames //
                        //
                        // TODO: add options to do intensity weighting or MPcorrection
                        // TODO: consider only doing the initialisation and data preparation once, currently repeated for every parameter --> speed improvement possible

                        eSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], nFramesOnGPU, nframeArray[nfi], blockSize, null, true, true, doFHTinterpolation, null, 0,1);
                        eSRRF.resetFramePosition();
                        eSRRF.loadDriftXYGPUbuffer(shiftXtempOdd, shiftYtempOdd);

                        userPressedEscape = calculateLiveSRRFframeLoad(imsAllRawData, nfi, 1, eSRRF);
                        if (userPressedEscape) {
                            eSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }

                        imsBuffer = eSRRF.imsSRRF;

                        fpOddAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                        fpOddVAR = imsBuffer.getProcessor(2).convertToFloatProcessor();
                        fpOddTAC2 = imsBuffer.getProcessor(3).convertToFloatProcessor();

                        // Calculate and get the reconstruction from the even frames // TODO: add options to do intensity weighting or MPcorrection
                        eSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], nFramesOnGPU, nframeArray[nfi], blockSize, null, true, true, doFHTinterpolation, null, 0, 1);
                        eSRRF.resetFramePosition();
                        eSRRF.loadDriftXYGPUbuffer(shiftXtempEven, shiftYtempEven);

                        calculateLiveSRRFframeLoad(imsAllRawData, nfi, 2, eSRRF);
                        imsBuffer = eSRRF.imsSRRF;

                        if (calculateAVG){
                            fpEvenAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                            pixelsFRCresolutionAVG[fwhmArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddAVG, fpEvenAVG, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        }
                        if (calculateVAR){
                            fpEvenVAR = imsBuffer.getProcessor(2).convertToFloatProcessor();
                            pixelsFRCresolutionVAR[fwhmArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddVAR, fpEvenVAR, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        }
                        if (calculateTAC2){
                            fpEvenTAC2 = imsBuffer.getProcessor(3).convertToFloatProcessor();
                            pixelsFRCresolutionTAC2[fwhmArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddTAC2, fpEvenTAC2, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        }

                        
                    } else {
                        // TODO: add options to do intensity weighting or MPcorrection
                        eSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], nFramesOnGPU, nframeArray[nfi], blockSize, null, true, true, doFHTinterpolation, null, 0, 1);
                        eSRRF.resetFramePosition();
                        eSRRF.loadDriftXYGPUbuffer(shiftXtemp, shiftYtemp);

                        calculateLiveSRRFframeLoad(imsAllRawData, nfi, 0, eSRRF);
                        imsBuffer = eSRRF.imsSRRF;
                    }

                    if (calculateAVG) imsSRRFavg.addSlice(label, imsBuffer.getProcessor(1));
                    if (calculateVAR) imsSRRFvar.addSlice(label, imsBuffer.getProcessor(2));
                    if (calculateTAC2) imsSRRFtac2.addSlice(label, imsBuffer.getProcessor(3));

                    // Interpolated image
                    imsInt.addSlice(label, imsBuffer.getProcessor(4));

                    // Error mapping
                    if (doErrorMapping) {

                        if (calculateAVG){
                            // Error maps for AVG
                            IJ.showStatus("Optimising Sigma...");
                            errorMapCalculator.optimise(imsBuffer.getProcessor(4), imsBuffer.getProcessor(1));

                            IJ.showStatus("Calculating error map...");
                            errorMapCalculator.calculateErrorMap();
                            pixelsRMSEavg[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                            pixelsPPMCCavg[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                            if (showErrorMaps) imsErrorMapAVG.addSlice(label, errorMapCalculator.fpErrorMap);
                            if (showRSC) imsRSCavg.addSlice(label, errorMapCalculator.fpSRC);
                        }

                        if (calculateVAR){
                            // Error maps for VAR
                            IJ.showStatus("Optimising Sigma...");
                            errorMapCalculator.optimise(imsBuffer.getProcessor(4), imsBuffer.getProcessor(2));

                            IJ.showStatus("Calculating error map...");
                            errorMapCalculator.calculateErrorMap();
                            pixelsRMSEvar[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                            pixelsPPMCCvar[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                            if (showErrorMaps) imsErrorMapVAR.addSlice(label, errorMapCalculator.fpErrorMap);
                            if (showRSC) imsRSCvar.addSlice(label, errorMapCalculator.fpSRC);
                        }

                        if (calculateTAC2){
                            // Error maps for TAC2
                            IJ.showStatus("Optimising Sigma...");
                            errorMapCalculator.optimise(imsBuffer.getProcessor(4), imsBuffer.getProcessor(3));

                            IJ.showStatus("Calculating error map...");
                            errorMapCalculator.calculateErrorMap();
                            pixelsRMSEtac2[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                            pixelsPPMCCtac2[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                            if (showErrorMaps) imsErrorMapTAC2.addSlice(label, errorMapCalculator.fpErrorMap);
                            if (showRSC) imsRSCtac2.addSlice(label, errorMapCalculator.fpSRC);
                        }
                    }


                    // Increment the reconstruction counter
                    r++;
                }
            }

            if (calculateRSE) {

                if (calculateAVG){
                    imsRMSEavg.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSEavg), nfi + 1);
                    imsRMSEavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateVAR){
                    imsRMSEvar.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSEvar), nfi + 1);
                    imsRMSEvar.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateTAC2){
                    imsRMSEtac2.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSEtac2), nfi + 1);
                    imsRMSEtac2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }

            if (calculateRSP) {

                if (calculateAVG){
                    imsPPMCCavg.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCCavg), nfi + 1);
                    imsPPMCCavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateVAR){
                    imsPPMCCvar.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCCvar), nfi + 1);
                    imsPPMCCvar.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateTAC2){
                    imsPPMCCtac2.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCCtac2), nfi + 1);
                    imsPPMCCtac2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }

            if (calculateFRC) {

                if (calculateAVG){
                    imsFRCresolutionAVG.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionAVG), nfi + 1);
                    imsFRCresolutionAVG.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateVAR){
                    imsFRCresolutionVAR.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionVAR), nfi + 1);
                    imsFRCresolutionVAR.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (calculateAVG){
                    imsFRCresolutionTAC2.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionTAC2), nfi + 1);
                    imsFRCresolutionTAC2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }
        }

        // Release the GPU
        eSRRF.release();

        //Display results
        Calibration sweepMapCalib = new Calibration();
        sweepMapCalib.setXUnit("Radius");
        sweepMapCalib.setYUnit("Sensitivity");

        if (fwhmArray.length > 1) {
            sweepMapCalib.pixelWidth = fwhmArray[1] - fwhmArray[0];
            sweepMapCalib.xOrigin = -fwhmArray[0] / (fwhmArray[1] - fwhmArray[0]);
        }
        else{
            sweepMapCalib.pixelWidth = fwhmArray[0];
            sweepMapCalib.xOrigin = -1;
        }

        if (sensitivityArray.length > 1) {
            sweepMapCalib.pixelHeight = sensitivityArray[1]-sensitivityArray[0];
            sweepMapCalib.yOrigin = -(float) sensitivityArray[0] / (float) (sensitivityArray[1]-sensitivityArray[0]);
        }
        else{
            sweepMapCalib.pixelHeight = sensitivityArray[0];
            sweepMapCalib.yOrigin = -1;
        }


        if (nframeArray.length > 1){
            sweepMapCalib.pixelDepth = nframeArray[1]-nframeArray[0];
            sweepMapCalib.zOrigin = -(float) nframeArray[0] / (float) (nframeArray[1]-nframeArray[0]);
        }
        else{
            sweepMapCalib.pixelDepth = nframeArray[0];
            sweepMapCalib.zOrigin = -1;
        }


        if (calculateAVG){
            displayImagePlus(imsSRRFavg, " - eSRRF (AVG)", cal, "");
//            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - eSRRF (AVG)", imsSRRFavg);
//            impSRRFavg.setCalibration(cal);
//            IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
//            impSRRFavg.show();
        }
        if (calculateVAR){
            displayImagePlus(imsSRRFvar, " - eSRRF (VAR)", cal, "");
//            ImagePlus impSRRFvar = new ImagePlus(imp.getTitle() + " - eSRRF (VAR)", imsSRRFvar);
//            impSRRFvar.setCalibration(cal);
//            IJ.run(impSRRFvar, "Enhance Contrast", "saturated=0.5");
//            impSRRFvar.show();
        }
        if (calculateTAC2){
            displayImagePlus(imsSRRFtac2, " - eSRRF (TAC2)", cal, "");
//            ImagePlus impSRRFtac2 = new ImagePlus(imp.getTitle() + " - eSRRF (TAC2)", imsSRRFtac2);
//            impSRRFtac2.setCalibration(cal);
//            IJ.run(impSRRFtac2, "Enhance Contrast", "saturated=0.5");
//            impSRRFtac2.show();
        }

        // Interpolated image
        displayImagePlus(imsInt, " - Interpolated image", cal, "");
//        ImagePlus impInt = new ImagePlus(imp.getTitle() + " - Interpolated image", imsInt);
//        impInt.setCalibration(cal);
//        IJ.run(impInt, "Enhance Contrast", "saturated=0.5");
//        impInt.show();

        // Showing error maps
        if (showErrorMaps) {

            if (calculateAVG){
                displayImagePlus(imsErrorMapAVG, " - ErrorMap (AVG)", cal, "ErrorMap-LUT");
//                ImagePlus impErrorMapAVG = new ImagePlus(imp.getTitle() + " - ErrorMap (AVG)", imsErrorMapAVG);
//                impErrorMapAVG.setCalibration(cal);
//                IJ.run(impErrorMapAVG, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impErrorMapAVG);
//                impErrorMapAVG.show();
            }
            if (calculateVAR){
                displayImagePlus(imsErrorMapVAR, " - ErrorMap (VAR)", cal, "ErrorMap-LUT");
//                ImagePlus impErrorMapVAR = new ImagePlus(imp.getTitle() + " - ErrorMap (VAR)", imsErrorMapVAR);
//                impErrorMapVAR.setCalibration(cal);
//                IJ.run(impErrorMapVAR, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impErrorMapVAR);
//                impErrorMapVAR.show();
            }
            if (calculateTAC2){
                displayImagePlus(imsErrorMapTAC2, " - ErrorMap (TAC2)", cal, "ErrorMap-LUT");
//                ImagePlus impErrorMapTAC2 = new ImagePlus(imp.getTitle() + " - ErrorMap (TAC2)", imsErrorMapTAC2);
//                impErrorMapTAC2.setCalibration(cal);
//                IJ.run(impErrorMapTAC2, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impErrorMapTAC2);
//                impErrorMapTAC2.show();
            }
        }

        if (showRSC) {

            if (calculateAVG){
                displayImagePlus(imsRSCavg, " - rescaled eSRRF (AVG)", cal, "");

//                ImagePlus impRSCavg = new ImagePlus(imp.getTitle() + " - rescaled eSRRF (AVG)", imsRSCavg);
//                impRSCavg.setCalibration(cal);
//                IJ.run(impRSCavg, "Enhance Contrast", "saturated=0.5");
//                impRSCavg.show();
            }
            if (calculateVAR){
                displayImagePlus(imsRSCvar, " - rescaled eSRRF (VAR)", cal, "");

//                ImagePlus impRSCvar = new ImagePlus(imp.getTitle() + " - rescaled eSRRF (VAR)", imsRSCvar);
//                impRSCvar.setCalibration(cal);
//                IJ.run(impRSCvar, "Enhance Contrast", "saturated=0.5");
//                impRSCvar.show();
            }
            if (calculateTAC2){
                displayImagePlus(imsRSCtac2, " - rescaled eSRRF (TAC2)", cal, "");

//                ImagePlus impRSCtac2 = new ImagePlus(imp.getTitle() + " - rescaled eSRRF (TAC2)", imsRSCtac2);
//                impRSCtac2.setCalibration(cal);
//                IJ.run(impRSCtac2, "Enhance Contrast", "saturated=0.5");
//                impRSCtac2.show();
            }
        }

        if (calculateRSE) {

            if (calculateAVG){
                displayImagePlus(imsRMSEavg, " - RMSE sweep map (AVG)", cal, "ErrorMap-LUT");

//                ImagePlus impRMSEavg = new ImagePlus(imp.getTitle() + " - RMSE sweep map (AVG)", imsRMSEavg);
//                impRMSEavg.setCalibration(sweepMapCalib.copy());
//                IJ.run(impRMSEavg, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impRMSEavg);
//                impRMSEavg.show();

                IJ.run("Maximize", "");
            }
            if (calculateVAR){
                displayImagePlus(imsRMSEvar, " - RMSE sweep map (VAR)", cal, "ErrorMap-LUT");

//                ImagePlus impRMSEvar = new ImagePlus(imp.getTitle() + " - RMSE sweep map (VAR)", imsRMSEvar);
//                impRMSEvar.setCalibration(sweepMapCalib.copy());
//                IJ.run(impRMSEvar, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impRMSEvar);
//                impRMSEvar.show();

                IJ.run("Maximize", "");
            }
            if (calculateTAC2){
                displayImagePlus(imsRMSEtac2, " - RMSE sweep map (TAC2)", cal, "ErrorMap-LUT");

//                ImagePlus impRMSEtac2 = new ImagePlus(imp.getTitle() + " - RMSE sweep map (TAC2)", imsRMSEtac2);
//                impRMSEtac2.setCalibration(sweepMapCalib.copy());
//                IJ.run(impRMSEtac2, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impRMSEtac2);
//                impRMSEtac2.show();

                IJ.run("Maximize", "");
            }
        }

        if (calculateRSP) {

            if (calculateAVG){
                displayImagePlus(imsPPMCCavg, " - RSP sweep map (AVG)", cal, "ErrorMap-LUT");

//                ImagePlus impPPMCCavg = new ImagePlus(imp.getTitle() + " - RSP sweep map (AVG)", imsPPMCCavg);
//                impPPMCCavg.setCalibration(sweepMapCalib.copy());
//                IJ.run(impPPMCCavg, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impPPMCCavg);
//                impPPMCCavg.show();

                IJ.run("Maximize", "");
            }
            if (calculateVAR){
                displayImagePlus(imsPPMCCvar, " - RSP sweep map (VAR)", cal, "ErrorMap-LUT");

//                ImagePlus impPPMCCvar = new ImagePlus(imp.getTitle() + " - RSP sweep map (VAR)", imsPPMCCvar);
//                impPPMCCvar.setCalibration(sweepMapCalib.copy());
//                IJ.run(impPPMCCvar, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impPPMCCvar);
//                impPPMCCvar.show();

                IJ.run("Maximize", "");
            }
            if (calculateTAC2){
                displayImagePlus(imsPPMCCtac2, " - RSP sweep map (TAC2)", cal, "ErrorMap-LUT");

//                ImagePlus impPPMCCtac2 = new ImagePlus(imp.getTitle() + " - RSP sweep map (TAC2)", imsPPMCCtac2);
//                impPPMCCtac2.setCalibration(sweepMapCalib.copy());
//                IJ.run(impPPMCCtac2, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impPPMCCtac2);
//                impPPMCCtac2.show();

                IJ.run("Maximize", "");
            }
        }

        if (calculateFRC) {

            if (calculateAVG){
                displayImagePlus(imsFRCresolutionAVG, " - FRC resolution sweep map (AVG)", cal, "FRC-LUT");

//                ImagePlus impFRCresolutionAVG = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (AVG)", imsFRCresolutionAVG);
//                impFRCresolutionAVG.setCalibration(sweepMapCalib.copy());
//                IJ.run(impFRCresolutionAVG, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_FRC(impFRCresolutionAVG);
//                impFRCresolutionAVG.show();

                IJ.run("Maximize", "");
            }
            if (calculateVAR){
                displayImagePlus(imsFRCresolutionVAR, " - FRC resolution sweep map (VAR)", cal, "FRC-LUT");

//                ImagePlus impFRCresolutionVAR = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (VAR)", imsFRCresolutionVAR);
//                impFRCresolutionVAR.setCalibration(sweepMapCalib.copy());
//                IJ.run(impFRCresolutionVAR, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_FRC(impFRCresolutionVAR);
//                impFRCresolutionVAR.show();

                IJ.run("Maximize", "");
            }
            if (calculateTAC2){
                displayImagePlus(imsFRCresolutionTAC2, " - FRC resolution sweep map (TAC2)", cal, "FRC-LUT");

//                ImagePlus impFRCresolutionTAC2 = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (TAC2)", imsFRCresolutionTAC2);
//                impFRCresolutionTAC2.setCalibration(sweepMapCalib.copy());
//                IJ.run(impFRCresolutionTAC2, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_FRC(impFRCresolutionTAC2);
//                impFRCresolutionTAC2.show();

                IJ.run("Maximize", "");
            }
        }

        // Run garbage collector
        System.gc();
        // Run garbage collector
        System.gc(); // twice is better than once !

        IJ.run("Cascade", "");

        IJ.log("-------------------------------------");
        IJ.log("RAM used: " + IJ.freeMemory());
        IJ.log("Bye-bye !");

    }

    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    //    --- Grab settings ---
    private boolean grabSettings(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();

        calculateAVG = gd.getNextBoolean();
        calculateVAR = gd.getNextBoolean();
        calculateTAC2 = gd.getNextBoolean();
        getInterpolatedImage = gd.getNextBoolean();

        correctVibration = gd.getNextBoolean();
        showRecons = gd.getNextBoolean();


        float fwhm0 = (float) gd.getNextNumber();
        float deltafwhm = (float) gd.getNextNumber();
        int n_fwhm = (int) gd.getNextNumber();

        int S0 = (int) gd.getNextNumber();
        int deltaS = (int) gd.getNextNumber();
        int n_S = (int) gd.getNextNumber();

        int nf0 = (int) gd.getNextNumber();
        int deltanf = (int) gd.getNextNumber();
        int n_nf = (int) gd.getNextNumber();

        // --- Prepare the variables ---
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
        if (calculateFRC) n_nfToUse = Math.min((nSlices / 2 - nf0) / deltanf + 1, n_nf);
        else n_nfToUse = Math.min((nSlices - nf0) / deltanf + 1, n_nf);

        if (n_nfToUse < 0){
            return false;
        }

        nframeArray = new int[n_nfToUse];
        for (int i = 0; i < n_nfToUse; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }

        boolean showAdvancedSettings = gd.getNextBoolean();
        if (showAdvancedSettings && !previousAdvSettings) advancedSettingsGUI();
        previousAdvSettings = showAdvancedSettings;

        prefs.set("magnification", magnification);
        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateVAR", calculateVAR);
        prefs.set("calculateTAC2", calculateTAC2);
        prefs.set("getInterpolatedImage", getInterpolatedImage);

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

        prefs.set("showRecons", showRecons);
        prefs.set("showErrorMaps", showErrorMaps);
        prefs.set("showRSC", showRSC);

        prefs.set("calculateRSE", calculateRSE);
        prefs.set("calculateRSP", calculateRSP);

        prefs.set("fixSigma", fixSigma);
        prefs.set("fixedSigma", fixedSigma);

        prefs.set("cropBorder", cropBorder);
        prefs.set("cropSize", cropSize);

        prefs.set("calculateFRC", calculateFRC);
        prefs.set("doFHTinterpolation", doFHTinterpolation);

        prefs.set("blockSize", blockSize);
        prefs.set("singleFrameLoad", singleFrameLoad);

        prefs.save();

        return true;
    }


    private boolean calculateLiveSRRFframeLoad(ImageStack imsAllRawData, int nf, int mode, LiveSRRF_CL eSRRF) {

//        ImageStack imsThisRawData;
        ImageStack imsThisRawData = new ImageStack(width, height);
        boolean userPressedEscape = false;
        int fmode;

        IJ.showStatus("Calculating eSRRF image...");
        for (int f = 1; f <= nframeArray[nf]; f++) {
            if (singleFrameLoad) imsThisRawData = new ImageStack(width, height); // re-initialise the stack every time
//            imsThisRawData = new ImageStack(width, height);

            if (mode == 0) fmode = f;  // no FRC
            else if (mode == 1) fmode = 2 * (f - 1) + 1; // FRC odd frames
            else fmode = 2 * f; // FRC even frames

            imsThisRawData.addSlice(imsAllRawData.getProcessor(fmode));
//            userPressedEscape = eSRRF.calculateSRRF(imsThisRawData);

            if (singleFrameLoad) {
                eSRRF.finishQueue(); // TODO: this seems necessary to provide reliable results, not sure why it's not necessary in the standard recon code though
                userPressedEscape = eSRRF.calculateSRRF(imsThisRawData);
                // Check if user is cancelling calculation
                if (userPressedEscape) {
                    return userPressedEscape;
                }
            }
            // Check if user is cancelling calculation
            if (userPressedEscape) {
                return userPressedEscape;
            }
        }


        // load the whole stack (hoping for the best in terms of memory)
        if (!singleFrameLoad) {
            userPressedEscape = eSRRF.calculateSRRF(imsThisRawData);
        }

        IJ.showStatus("Reading eSRRF image...");
        eSRRF.readSRRFbuffer();

        return userPressedEscape;

    }

    private void advancedSettingsGUI(){

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        GenericDialog gd = new GenericDialog("eSRRF Parameter sweep - Advanced settings");

        gd.addMessage("-=-= Error maps =-=-\n", headerFont);
        gd.addCheckbox("Show all error maps (default: on)", prefs.get("showErrorMaps", true));
        gd.addToSameRow();
        gd.addCheckbox("Show rescaled reconstructions (default: off)", prefs.get("showRSC", false));

        gd.addCheckbox("Calculate RSE (default: off)", prefs.get("calculateRSE", false));
        gd.addToSameRow();
        gd.addCheckbox("Calculate RSP (default: off)", prefs.get("calculateRSP", false));

        gd.addCheckbox("Fix sigma (default: off)", prefs.get("fixSigma", false));
        gd.addToSameRow();
        gd.addNumericField("Sigma (in pixels, used if fixed)", prefs.get("fixedSigma", 1), 2);
        gd.addMessage("Sigma can be evaluated by: (0.21 x Emission wavelength (in nm) / NA) / pixel size (in nm)\n"+
                "It should typically be around ~1 pixel.");

        gd.addCheckbox("Crop border (default: on)", prefs.get("cropBorder", true));
        gd.addToSameRow();
        gd.addNumericField("Pixels to crop (in pixels, used if selected)", prefs.get("cropSize", 3), 0);
        gd.addMessage("Cropping the border of the images avoids biasing error mapping analysis if extrapolation\n"+
                " artefacts may occur.");

        gd.addMessage("-=-= FRC resolution =-=-\n", headerFont);
        gd.addCheckbox("Calculate FRC (default: off)", prefs.get("calculateFRC", false));
        gd.addMessage("Calculating FRC will split all dataset in two halves and therefore the maximum number\n" +
                "of frames will be half of the total frames available in the dataset.");

        gd.addMessage("-=-= Reconstruction settings =-=-\n", headerFont);
        gd.addCheckbox("Use FHT interpolation", prefs.get("doFHTinterpolation",true));

        gd.addMessage("-=-= GPU processing =-=-\n", headerFont);
        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use more resources and\n" +
                " may slow down your computer.");
        gd.addCheckbox("Enable single frame loading", prefs.get("singleFrameLoad",true));
        gd.addMessage("Single frame loading is slower but is more stable due to lower GPU memory requirements.");


        gd.showDialog();

        if (!gd.wasCanceled()){
            showErrorMaps = gd.getNextBoolean();
            showRSC = gd.getNextBoolean();

            calculateRSE = gd.getNextBoolean();
            calculateRSP = gd.getNextBoolean();

            // If any of these things is ticked then do the analysis
            doErrorMapping = (showErrorMaps || showRSC || calculateRSE || calculateRSP);

            fixSigma = gd.getNextBoolean();
            fixedSigma = (float) gd.getNextNumber();

            cropBorder = gd.getNextBoolean();
            cropSize = (int) gd.getNextNumber();

            calculateFRC = gd.getNextBoolean();
            doFHTinterpolation = gd.getNextBoolean();

            blockSize = (int) gd.getNextNumber();
            singleFrameLoad = gd.getNextBoolean();
        }
    }

    // --- Main GUI Dialog listener ---
    class MyDialogListenerMainGUI implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {
            return grabSettings(gd);
        }
    }

    // ---- Displayer!! ----
    private void displayImagePlus(ImageStack ims, String titleAppendix, Calibration cal, String nameLUT) {

        ImagePlus imp = new ImagePlus(imageTitle + titleAppendix, ims);
        imp.setCalibration(cal.copy());
        IJ.run(imp, "Enhance Contrast", "saturated=0.5");
        if (nameLUT.equals("ErrorMap-LUT")) applyLUT_SQUIRREL_Errors(imp);
        else if (nameLUT.equals("FRC-LUT")) applyLUT_SQUIRREL_FRC(imp);

        imp.show();
    }


}
