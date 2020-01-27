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
            nRecons,
            nSlices,
            width,
            height,
            blockSize,
            cropSize;

    private boolean showRecons,
    //            calculateAVG,
//            calculateVAR,
//            calculateTAC2,
    showErrorMaps,
            showRSC,
            calculateFRC,
            calculateRSE,
            calculateRSP,
            doErrorMapping,
            fixSigma,
            correctVibration,
            cropBorder;

    private String[] reconsNames;

    private boolean[] calculateReconArray;

    private float[] radiusArray;
    private int[] sensitivityArray,
            nframeArray;

    private float fixedSigma;

    private final String LiveSRRFVersion = "v1.2d-fhi.3";
    private float[] shiftX, shiftY;

    private String imageTitle;

//    private String chosenTemporalAnalysis;

    // Advanced formats
    private final NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
//    private liveSRRF_CL liveSRRF;

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
        IJ.log("LiveSRRF - Parameters sweep " + LiveSRRFVersion);
        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));

        // Initialise the liveSRRF engine (creates a context really)
        LiveSRRF_CL liveSRRF = new LiveSRRF_CL();
        nRecons = liveSRRF.nReconstructions;
        reconsNames = liveSRRF.reconNames;
        calculateReconArray = new boolean[nRecons];

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("LiveSRRF - Parameters sweep " + LiveSRRFVersion);
        gd.addMessage("-=-= LiveSRRF reconstruction =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 4)", prefs.get("magnification", 4), 0);

        for (int i = 0; i < nRecons; i++) {
            gd.addCheckbox(reconsNames[i]+" reconstruction (default: on)", prefs.get("calculate"+reconsNames[i], true));
        }
//        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
//        gd.addCheckbox("VAR reconstruction (default: off)", prefs.get("calculateVAR", false));
//        gd.addCheckbox("TAC2 reconstruction (default: off)", prefs.get("calculateTAC2", false));

        gd.addCheckbox("Correct vibration (default: off)", prefs.get("correctVibration", false));
        gd.addCheckbox("Show all reconstructions (default: on)", prefs.get("showRecons", true));

        gd.addMessage("-=-= Sweeping LiveSRRF parameters =-=-\n", headerFont);
        gd.addMessage("Radius\n", headerFont);
        gd.addNumericField("Start", prefs.get("radius0", 2), 2);
        gd.addToSameRow();
        gd.addNumericField("Delta", prefs.get("deltaRadius", 0.5f), 2);
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

        gd.addMessage("-=-= GPU processing =-=-\n", headerFont);
        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use more resources and\n" +
                " may slow down your computer.");

        gd.addHelp("https://www.youtube.com/watch?v=Vs6awg-BJHo"); // If you don't know how to sweep
        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            liveSRRF.release();
            return;
        }

        boolean goodToGo = grabSettings(gd);
        if (!goodToGo){
            liveSRRF.release();
            IJ.log("-------------------------------------");
            IJ.log("Aborted: Not enough frames in the stack for this analysis.");
            return;
        }


        // Log the settings used
        IJ.log("Magnification: " + magnification);
        if (correctVibration) IJ.log("Vibration correction: on");
        else IJ.log("Vibration correction: off");

        if (doErrorMapping) {
            IJ.log("Error mapping parameters:");
            if (fixSigma) IJ.log("Sigma is fixed to " + fixedSigma + " pixels.");
            else IJ.log("Sigma is optimised for each reconstructions.");
            if (cropBorder) IJ.log("The borders will be cropped by "+cropSize+" pixels");
        }

        int n_calculation = nframeArray.length * sensitivityArray.length * radiusArray.length;
        IJ.log("Number of calculations planned: " + n_calculation);

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

        ImageStack[] imsSRRFarray = new ImageStack[nRecons];
        ImageStack[] imsErrorMapArray = new ImageStack[nRecons];
        ImageStack[] imsRSCarray = new ImageStack[nRecons];


        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);

        for (int i = 0; i < nRecons; i++) {
            imsSRRFarray[i] = new ImageStack(width * magnification, height * magnification);
        }
//        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
//        ImageStack imsSRRFvar = new ImageStack(width * magnification, height * magnification);
//        ImageStack imsSRRFtac2 = new ImageStack(width * magnification, height * magnification);

        // Initialising the Error maps variables
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF(imp, magnification, fixSigma, fixedSigma, cropBorder, cropSize);

        float[][] pixelsRMSEarray = null;
        float[][] pixelsPPMCarray = null;
        float[][] pixelsFRCresolutionArray = null;

//        float[] pixelsRMSEavg = null;
//        float[] pixelsPPMCCavg = null;
//        float[] pixelsRMSEvar = null;
//        float[] pixelsPPMCCvar = null;
//        float[] pixelsRMSEtac2 = null;
//        float[] pixelsPPMCCtac2 = null;

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

//        ImageStack imsErrorMapAVG = new ImageStack(widthSquirrel, heightSquirrel);
//        ImageStack imsErrorMapSTD = new ImageStack(widthSquirrel, heightSquirrel);
//        ImageStack imsErrorMapTAC2 = new ImageStack(widthSquirrel, heightSquirrel);

//        ImageStack imsRSCavg = new ImageStack(widthSquirrel, heightSquirrel);
//        ImageStack imsRSCvar = new ImageStack(widthSquirrel, heightSquirrel);
//        ImageStack imsRSCtac2 = new ImageStack(widthSquirrel, heightSquirrel);


        ImageStack[] imsRMSEarray = new ImageStack[nRecons];
        ImageStack[] imsPPMCCarray = new ImageStack[nRecons];

        // These images are small so it doesn't matter if they are initialised in any case
        for (int i = 0; i < nRecons; i++) {
            imsErrorMapArray[i] = new ImageStack(widthSquirrel, heightSquirrel);
            imsRSCarray[i] = new ImageStack(widthSquirrel, heightSquirrel);
            imsRMSEarray[i] = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
            imsPPMCCarray[i] = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
        }
//        ImageStack imsRMSEavg = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsPPMCCavg = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsRMSEvar = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsPPMCCvar = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsRMSEtac2 = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsPPMCCtac2 = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);

        // Initialising the shift correction variables
        int maxnFrame;
        if (calculateFRC) maxnFrame = 2 * nframeArray[nframeArray.length - 1];
        else maxnFrame = nframeArray[nframeArray.length - 1];

        shiftX = new float[maxnFrame];
        shiftY = new float[maxnFrame];
        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp);

        if (correctVibration) {
            shiftCalculator.calculateShiftArray(1, maxnFrame);
            shiftX = shiftCalculator.shiftX;
            shiftY = shiftCalculator.shiftY;
        }

        float[][] shiftXYtemp;
//        float[] shiftYtemp;

        float[][] shiftXYtempOdd;
//        float[] shiftYtempOdd;

        float[][] shiftXYtempEven;
//        float[] shiftYtempEven;


        // Initialising the FRC variables
        FRC frcCalculator = new FRC();

        FloatProcessor[][] fpOddEvenArray = new FloatProcessor[nRecons][2];
//        FloatProcessor fpOddAVG;
//        FloatProcessor fpEvenAVG;
//        FloatProcessor fpOddVAR;
//        FloatProcessor fpEvenVAR;
//        FloatProcessor fpOddTAC2;
//        FloatProcessor fpEvenTAC2;

//        float[] pixelsFRCresolutionAVG = null;
//        float[] pixelsFRCresolutionVAR = null;
//        float[] pixelsFRCresolutionTAC2 = null;

        // These images are small so it doesn't matter if they are initialised in any case
        ImageStack[] imsFRCresolutionArray = new ImageStack[nRecons];
        for (int i = 0; i < nRecons; i++) {
            imsFRCresolutionArray[i] = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
        }

//        ImageStack imsFRCresolutionAVG = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsFRCresolutionVAR = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);
//        ImageStack imsFRCresolutionTAC2 = new ImageStack(radiusArray.length, sensitivityArray.length, nframeArray.length);

        // Set the number of reconstructions
        int r = 0;
        boolean userPressedEscape;

        for (int nfi = 0; nfi < nframeArray.length; nfi++) {

            if (doErrorMapping){
                pixelsRMSEarray = new float[nRecons][(radiusArray.length) * (sensitivityArray.length)];
                pixelsPPMCarray = new float[nRecons][(radiusArray.length) * (sensitivityArray.length)];
            }
//            if (doErrorMapping) {
//                pixelsRMSEavg = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsPPMCCavg = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsRMSEvar = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsPPMCCvar = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsRMSEtac2 = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsPPMCCtac2 = new float[(radiusArray.length) * (sensitivityArray.length)];
//            }
//
            if (calculateFRC) pixelsFRCresolutionArray = new float[nRecons][(radiusArray.length) * (sensitivityArray.length)];

//            if (calculateFRC) {
//                pixelsFRCresolutionAVG = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsFRCresolutionVAR = new float[(radiusArray.length) * (sensitivityArray.length)];
//                pixelsFRCresolutionTAC2 = new float[(radiusArray.length) * (sensitivityArray.length)];
//            }


            shiftXYtemp = new float[nframeArray[nfi]][2];
//            shiftYtemp = new float[nframeArray[nfi]];

            shiftXYtempOdd = new float[nframeArray[nfi]][2];
//            shiftYtempOdd = new float[nframeArray[nfi]];
            shiftXYtempEven = new float[nframeArray[nfi]][2];
//            shiftYtempEven = new float[nframeArray[nfi]];

            // Re-adjust the array of shifts depending on the frames used
            if (calculateFRC) {
                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXYtempOdd[i][0] = shiftX[2 * i + 1];
                    shiftXYtempOdd[i][1] = shiftY[2 * i + 1];
//                IJ.log("ShiftX="+shiftXtemp[i]);
//                    shiftYtempOdd[i] = shiftY[2 * i + 1];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }

                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXYtempEven[i][0] = shiftX[2 * i];
                    shiftXYtempEven[i][1] = shiftY[2 * i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
//                    shiftYtempEven[i] = shiftY[2 * i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }

            } else {
                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXYtemp[i][0] = shiftX[i];
                    shiftXYtemp[i][1] = shiftY[i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
//                    shiftYtemp[i] = shiftY[i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }
            }


            for (int si = 0; si < sensitivityArray.length; si++) {
                for (int fi = 0; fi < radiusArray.length; fi++) {

                    IJ.log("--------");
                    IJ.log("LiveSRRF frame: " + (r + 1) + "/" + n_calculation);
                    IJ.showProgress(r, n_calculation);

                    // Check if user is cancelling calculation
                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        liveSRRF.release();
                        IJ.log("-------------------------------------");
                        IJ.log("Reconstruction aborted by user.");
                        return;
                    }

                    IJ.log("Number of frame for SRRF: " + nframeArray[nfi]);
                    IJ.log("Radius: " + radiusArray[fi] + " pixels");
                    IJ.log("Sensitivity: " + sensitivityArray[si]);

                    String label = "R=" + radiusArray[fi] + "/S=" + sensitivityArray[si] + "/#fr=" + nframeArray[nfi];


                    // FRC resolution estimation
                    if (calculateFRC) {
                        // Calculate and get the reconstruction from the odd frames // TODO: add options to do intensity weighting or MPcorrection?
                        liveSRRF.initialise(width, height, magnification, radiusArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, false, "RobX");
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXYtempOdd);

                        userPressedEscape = calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 1, liveSRRF);
                        if (userPressedEscape) {
                            liveSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }

                        imsBuffer = liveSRRF.imsSRRF;

//                        fpOddAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
//                        fpOddVAR = imsBuffer.getProcessor(2).convertToFloatProcessor();
//                        fpOddTAC2 = imsBuffer.getProcessor(2).convertToFloatProcessor();
                        for (int i = 0; i < nRecons; i++) {
                            if (calculateReconArray[i]) fpOddEvenArray[i][0] = imsBuffer.getProcessor(i+1).convertToFloatProcessor();
                        }

                        // Calculate and get the reconstruction from the even frames // TODO: add options to do intensity weighting or MPcorrection
                        liveSRRF.initialise(width, height, magnification, radiusArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, true, "RobX");
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXYtempEven);

                        userPressedEscape = calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 2, liveSRRF);
                        if (userPressedEscape) {
                            liveSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }

                        imsBuffer = liveSRRF.imsSRRF;

                        for (int i = 0; i < nRecons; i++) {
                            if (calculateReconArray[i]) {
                                fpOddEvenArray[i][1] = imsBuffer.getProcessor(i+1).convertToFloatProcessor();
                                pixelsFRCresolutionArray[i][radiusArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddEvenArray[i][0], fpOddEvenArray[i][1], FRC.ThresholdMethod.FIXED_1_OVER_7);
                            }
                        }
////                        if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                        if (calculateReconArray[0]) {
//                            fpEvenAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
//                            pixelsFRCresolutionAVG[radiusArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddAVG, fpEvenAVG, FRC.ThresholdMethod.FIXED_1_OVER_7);
//                        }
////                        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                        if (calculateReconArray[1]) {
//                            fpEvenVAR = imsBuffer.getProcessor(2).convertToFloatProcessor();
//                            pixelsFRCresolutionVAR[radiusArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddVAR, fpEvenVAR, FRC.ThresholdMethod.FIXED_1_OVER_7);
//                        }
//
//                        if (calculateReconArray[2]) {
//                            fpEvenTAC2 = imsBuffer.getProcessor(3).convertToFloatProcessor();
//                            pixelsFRCresolutionTAC2[radiusArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddTAC2, fpEvenTAC2, FRC.ThresholdMethod.FIXED_1_OVER_7);
//                        }
                    }


                    else { // if (calculateFRC)
                        // TODO: add options to do intensity weighting or MPcorrection?
                        liveSRRF.initialise(width, height, magnification, radiusArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, true, "RobX");
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXYtemp);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 0, liveSRRF);
                        imsBuffer = liveSRRF.imsSRRF;
                    }


                    // Collate the reconstructions
//                    if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2]))
//                    if (calculateReconArray[0]) imsSRRFavg.addSlice(label, imsBuffer.getProcessor(1));
////                    if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2]))
//                    if (calculateReconArray[1]) imsSRRFvar.addSlice(label, imsBuffer.getProcessor(2));
//                    if (calculateReconArray[2]) imsSRRFtac2.addSlice(label, imsBuffer.getProcessor(3));

                    for (int i = 0; i < nRecons; i++) {
                        if (calculateReconArray[i]) imsSRRFarray[i].addSlice(label, imsBuffer.getProcessor(i+1));
                    }

                    // Interpolated image
                    imsInt.addSlice(label, imsBuffer.getProcessor(nRecons+1));

                    // Error mapping
                    if (doErrorMapping) {

                        for (int i = 0; i < nRecons; i++) {
                            if (calculateReconArray[i]){
                                // Error maps for AVG
                                IJ.showStatus("Optimising Sigma...");
                                errorMapCalculator.optimise(imsBuffer.getProcessor(nRecons+1), imsBuffer.getProcessor(i+1));

                                IJ.showStatus("Calculating error map...");
                                errorMapCalculator.calculateErrorMap();
                                pixelsRMSEarray[i][radiusArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                                pixelsPPMCarray[i][radiusArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                                if (showErrorMaps) imsErrorMapArray[i].addSlice(label, errorMapCalculator.fpErrorMap);
                                if (showRSC) imsRSCarray[i].addSlice(label, errorMapCalculator.fpSRC);
                            }
                        }

////                        if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                        if (calculateReconArray[0]){
//                            // Error maps for AVG
//                            IJ.showStatus("Optimising Sigma...");
//                            errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(1));
//
//                            IJ.showStatus("Calculating error map...");
//                            errorMapCalculator.calculateErrorMap();
//                            pixelsRMSEavg[radiusArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
//                            pixelsPPMCCavg[radiusArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;
//
//                            if (showErrorMaps) imsErrorMapAVG.addSlice(label, errorMapCalculator.fpErrorMap);
//                            if (showRSC) imsRSCavg.addSlice(label, errorMapCalculator.fpSRC);
//
//                        }
//
////                        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                        if (calculateReconArray[1]){
//
//                            // Error maps for STD
//                            IJ.showStatus("Optimising Sigma...");
//                            errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(2));
//
//                            IJ.showStatus("Calculating error map...");
//                            errorMapCalculator.calculateErrorMap();
//                            pixelsRMSEvar[radiusArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
//                            pixelsPPMCCvar[radiusArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;
//
//                            if (showErrorMaps) imsErrorMapSTD.addSlice(label, errorMapCalculator.fpErrorMap);
//                            if (showRSC) imsRSCvar.addSlice(label, errorMapCalculator.fpSRC);
//                        }
//                        //                        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                        if (calculateReconArray[2]){
//
//                            // Error maps for STD
//                            IJ.showStatus("Optimising Sigma...");
//                            errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(2));
//
//                            IJ.showStatus("Calculating error map...");
//                            errorMapCalculator.calculateErrorMap();
//                            pixelsRMSEtac2[radiusArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
//                            pixelsPPMCCtac2[radiusArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;
//
//                            if (showErrorMaps) imsErrorMapTAC2.addSlice(label, errorMapCalculator.fpErrorMap);
//                            if (showRSC) imsRSCtac2.addSlice(label, errorMapCalculator.fpSRC);
//                        }
                    }


                    // Increment the reconstruction counter
                    r++;
                }
            }

            if (calculateRSE) {

                for (int i = 0; i < nRecons; i++) {
                    if (calculateReconArray[i]){
                        imsRMSEarray[i].setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsRMSEarray[i]), nfi + 1);
                        imsRMSEarray[i].setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                    }
                }
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateAVG){
//                    imsRMSEavg.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsRMSEavg), nfi + 1);
//                    imsRMSEavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
//
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateVAR){
//                    imsRMSEvar.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsRMSEvar), nfi + 1);
//                    imsRMSEvar.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
//
//                //                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateTAC2){
//                    imsRMSEtac2.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsRMSEtac2), nfi + 1);
//                    imsRMSEtac2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
            }

            if (calculateRSP) {

                for (int i = 0; i < nRecons; i++) {
                    if (calculateReconArray[i]){
                        imsPPMCCarray[i].setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsPPMCarray[i]), nfi + 1);
                        imsPPMCCarray[i].setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                    }
                }
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateAVG){
//                    imsPPMCCavg.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsPPMCCavg), nfi + 1);
//                    imsPPMCCavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateVAR){
//                    imsPPMCCvar.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsPPMCCvar), nfi + 1);
//                    imsPPMCCvar.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
//
//                //                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateVAR){
//                    imsPPMCCtac2.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsPPMCCtac2), nfi + 1);
//                    imsPPMCCtac2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
            }

            if (calculateFRC) {
                for (int i = 0; i < nRecons; i++) {
                    if (calculateReconArray[i]){
                        imsFRCresolutionArray[i].setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsFRCresolutionArray[i]), nfi + 1);
                        imsFRCresolutionArray[i].setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                    }
                }
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateAVG){
//                    imsFRCresolutionAVG.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsFRCresolutionAVG), nfi + 1);
//                    imsFRCresolutionAVG.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
////                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateVAR){
//                    imsFRCresolutionVAR.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsFRCresolutionVAR), nfi + 1);
//                    imsFRCresolutionVAR.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
//
//                //                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//                if (calculateTAC2){
//                    imsFRCresolutionTAC2.setProcessor(new FloatProcessor(radiusArray.length, sensitivityArray.length, pixelsFRCresolutionTAC2), nfi + 1);
//                    imsFRCresolutionTAC2.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
//                }
            }
        }

        // Release the GPU
        liveSRRF.release();

        //Display results
        Calibration sweepMapCalib = new Calibration();
        sweepMapCalib.setXUnit("Radius");
        sweepMapCalib.setYUnit("Sensitivity");

        if (radiusArray.length > 1) {
            sweepMapCalib.pixelWidth = radiusArray[1] - radiusArray[0];
            sweepMapCalib.xOrigin = -radiusArray[0] / (radiusArray[1] - radiusArray[0]);
        }
        else{
            sweepMapCalib.pixelWidth = radiusArray[0];
            sweepMapCalib.xOrigin = -1;
        }

        if (sensitivityArray.length > 1) {
            sweepMapCalib.pixelHeight = sensitivityArray[1]-sensitivityArray[0];
            sweepMapCalib.yOrigin = -(float)sensitivityArray[0]/ (float)(sensitivityArray[1]-sensitivityArray[0]);
        }
        else{
            sweepMapCalib.pixelHeight = sensitivityArray[0];
            sweepMapCalib.yOrigin = -1;
        }


        if (nframeArray.length > 1){
            sweepMapCalib.pixelDepth = nframeArray[1]-nframeArray[0];
            sweepMapCalib.zOrigin = -(float)nframeArray[0]/ (float)(nframeArray[1]-nframeArray[0]);
        }
        else{
            sweepMapCalib.pixelDepth = nframeArray[0];
            sweepMapCalib.zOrigin = -1;
        }


        for (int i = 0; i < nRecons; i++) {
            if (calculateReconArray[i]) displayImagePlus(imsSRRFarray[i], " - LiveSRRF ("+reconsNames[i]+")", cal, "");
        }

//        // liveSRRF (AVG) reconstruction
//        if (calculateAVG){
//            displayImagePlus(imsSRRFavg, " - LiveSRRF (AVG)", cal, "");
////            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - LiveSRRF (AVG)", imsSRRFavg);
////            impSRRFavg.setCalibration(cal);
////            IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
////            impSRRFavg.show();
//        }
//
//        // liveSRRF (VAR) reconstruction
//        if (calculateVAR){
//            displayImagePlus(imsSRRFvar, " - LiveSRRF (VAR)", cal, "");
//
////            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - LiveSRRF (VAR)", imsSRRFvar);
////            impSRRFstd.setCalibration(cal);
////            IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
////            impSRRFstd.show();
//        }
//
//        // liveSRRF (VAR) reconstruction
//        if (calculateTAC2){
//            displayImagePlus(imsSRRFtac2, " - LiveSRRF (TAC2)", cal, "");
//
////            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - LiveSRRF (VAR)", imsSRRFvar);
////            impSRRFstd.setCalibration(cal);
////            IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
////            impSRRFstd.show();
//        }

//        // Interpolated image
//        ImagePlus impInt = new ImagePlus(imp.getTitle() + " - Interpolated image", imsInt);
//        impInt.setCalibration(cal);
//        IJ.run(impInt, "Enhance Contrast", "saturated=0.5");
//        impInt.show();

        displayImagePlus(imsInt, " - Interpolated image", cal, "");

        // Showing error maps
        if (showErrorMaps) {

            for (int i = 0; i < nRecons; i++) {
                if (calculateReconArray[i]) displayImagePlus(imsErrorMapArray[i], " - Error map ("+reconsNames[i]+")", cal, "ErrorMap-LUT");
            }

//            if (calculateReconArray[0]){
//                ImagePlus impErrorMapAVG = new ImagePlus(imp.getTitle() + " - ErrorMap (AVG)", imsErrorMapAVG);
//                impErrorMapAVG.setCalibration(cal);
//                IJ.run(impErrorMapAVG, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impErrorMapAVG);
//                impErrorMapAVG.show();
//            }
//            if (calculateReconArray[1]){
//                ImagePlus impErrorMapVAR = new ImagePlus(imp.getTitle() + " - ErrorMap (VAR)", imsErrorMapSTD);
//                impErrorMapVAR.setCalibration(cal);
//                IJ.run(impErrorMapVAR, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impErrorMapVAR);
//                impErrorMapVAR.show();
//            }
        }

        if (showRSC) {

            for (int i = 0; i < nRecons; i++) {
                if (calculateReconArray[i]) displayImagePlus(imsRSCarray[i], " - rescaled LiveSRRF ("+reconsNames[i]+")", cal, "");
            }
//            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateAVG){
//                ImagePlus impRSCavg = new ImagePlus(imp.getTitle() + " - rescaled liveSRRF (AVG)", imsRSCavg);
//                impRSCavg.setCalibration(cal);
//                IJ.run(impRSCavg, "Enhance Contrast", "saturated=0.5");
//                impRSCavg.show();
//            }
////            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateVAR){
//                ImagePlus impRSCstd = new ImagePlus(imp.getTitle() + " - rescaled liveSRRF (STD)", imsRSCvar);
//                impRSCstd.setCalibration(cal);
//                IJ.run(impRSCstd, "Enhance Contrast", "saturated=0.5");
//                impRSCstd.show();
//            }
        }

        if (calculateRSE) {
            for (int i = 0; i < nRecons; i++) {
                if (calculateReconArray[i]) displayImagePlus(imsRMSEarray[i], " - RMSE sweep map ("+reconsNames[i]+")", sweepMapCalib, "ErrorMap-LUT");
            }

////            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateAVG){
//
//                ImagePlus impRMSEavg = new ImagePlus(imp.getTitle() + " - RMSE sweep map (AVG)", imsRMSEavg);
//                impRMSEavg.setCalibration(sweepMapCalib.copy());
//                IJ.run(impRMSEavg, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impRMSEavg);
//                impRMSEavg.show();
//            }
////            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateVAR){
//
//                ImagePlus impRMSEstd = new ImagePlus(imp.getTitle() + " - RMSE sweep map (STD)", imsRMSEvar);
//                impRMSEstd.setCalibration(sweepMapCalib.copy());
//                IJ.run(impRMSEstd, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impRMSEstd);
//                impRMSEstd.show();
//            }
        }

        if (calculateRSP) {

            for (int i = 0; i < nRecons; i++) {
                if (calculateReconArray[i]) displayImagePlus(imsPPMCCarray[i], " - RSP sweep map ("+reconsNames[i]+")", sweepMapCalib, "ErrorMap-LUT");
            }

//            if (calculateAVG){
//
//                ImagePlus impPPMCCavg = new ImagePlus(imp.getTitle() + " - RSP sweep map (AVG)", imsPPMCCavg);
//                impPPMCCavg.setCalibration(sweepMapCalib.copy());
//                IJ.run(impPPMCCavg, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impPPMCCavg);
//                impPPMCCavg.show();
//            }
////            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateVAR){
//
//                ImagePlus impPPMCCstd = new ImagePlus(imp.getTitle() + " - RSP sweep map (VAR)", imsPPMCCvar);
//                impPPMCCstd.setCalibration(sweepMapCalib.copy());
//                IJ.run(impPPMCCstd, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_Errors(impPPMCCstd);
//                impPPMCCstd.show();
//            }
        }

        if (calculateFRC) {

            for (int i = 0; i < nRecons; i++) {
                if (calculateReconArray[i]) displayImagePlus(imsFRCresolutionArray[i], " - FRC resolution sweep map ("+reconsNames[i]+")", sweepMapCalib, "FRC-LUT");
            }
////            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateAVG){
//
//                ImagePlus impFRCresolutionAVG = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (AVG)", imsFRCresolutionAVG);
//                impFRCresolutionAVG.setCalibration(sweepMapCalib.copy());
//                IJ.run(impFRCresolutionAVG, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_FRC(impFRCresolutionAVG);
//                impFRCresolutionAVG.show();
//            }
////            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
//            if (calculateVAR){
//
//                ImagePlus impFRCresolutionSTD = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (VAR)", imsFRCresolutionVAR);
//                impFRCresolutionSTD.setCalibration(sweepMapCalib.copy());
//                IJ.run(impFRCresolutionSTD, "Enhance Contrast", "saturated=0.5");
//                applyLUT_SQUIRREL_FRC(impFRCresolutionSTD);
//                impFRCresolutionSTD.show();
//            }
        }

        // Run garbage collector
        System.gc();
        // Run garbage collector
        System.gc(); // twice is better than once !

        IJ.log("-------------------------------------");
        IJ.log("RAM used: " + IJ.freeMemory());
        IJ.log("Bye-bye !");



    }

//    -------------------------------------------------------------------------------------
//    -------------------------- Here lie dragons and functions ---------------------------
//    -------------------------------------------------------------------------------------


    //    --- Grab settings ---
    private boolean grabSettings(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();
//        chosenTemporalAnalysis = gd.getNextChoice();

        for (int i = 0; i < nRecons; i++) {
            calculateReconArray[i] = gd.getNextBoolean();
        }
//        calculateAVG = gd.getNextBoolean();
//        calculateVAR = gd.getNextBoolean();
//        calculateTAC2 = gd.getNextBoolean();

        correctVibration = gd.getNextBoolean();
        showRecons = gd.getNextBoolean();


        float radius0 = (float) gd.getNextNumber();
        float deltaRadius = (float) gd.getNextNumber();
        int n_fwhm = (int) gd.getNextNumber();

        int S0 = (int) gd.getNextNumber();
        int deltaS = (int) gd.getNextNumber();
        int n_S = (int) gd.getNextNumber();

        int nf0 = (int) gd.getNextNumber();
        int deltanf = (int) gd.getNextNumber();
        int n_nf = (int) gd.getNextNumber();


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

        blockSize = (int) gd.getNextNumber();

        radiusArray = new float[n_fwhm];
        for (int i = 0; i < n_fwhm; i++) {
            radiusArray[i] = radius0 + i * deltaRadius;
        }

        sensitivityArray = new int[n_S];
        for (int i = 0; i < n_S; i++) {
            sensitivityArray[i] = S0 + i * deltaS;
        }

        // Check that nf does not exceed nSlices
        int n_nfToUse;
        if (calculateFRC) n_nfToUse = Math.min((nSlices / 2 - nf0) / deltanf + 1, n_nf);
        else n_nfToUse = Math.min((nSlices - nf0) / deltanf + 1, n_nf);
        // TODO: ABORT IF NOT enough frames in original stack for even 1 analysis

        if (n_nfToUse < 0){
            return false;
        }

        nframeArray = new int[n_nfToUse];
        for (int i = 0; i < n_nfToUse; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }

        prefs.set("magnification", magnification);
//        prefs.set("chosenTemporalAnalysis", chosenTemporalAnalysis);

        for (int i = 0; i < nRecons; i++) {
            prefs.set("calculate"+reconsNames[i], calculateReconArray[i]);
        }
//        prefs.set("calculateAVG", calculateAVG);
//        prefs.set("calculateVAR", calculateVAR);
//        prefs.set("calculateTAC2", calculateTAC2);

        prefs.set("correctVibration", correctVibration);

        prefs.set("radius0", radius0);
        prefs.set("deltaRadius", deltaRadius);
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

        prefs.set("blockSize", blockSize);


        prefs.save();

        return true;
    }


    private boolean calculateLiveSRRFsingleframeLoad(ImageStack imsAllRawData, int nf, int mode, LiveSRRF_CL liveSRRF) {
        ImageStack imsThisRawData;
        boolean userPressedEscape = false;
        int fmode;

        IJ.showStatus("Calculating LiveSRRF image...");
        for (int f = 1; f <= nframeArray[nf]; f++) {
            imsThisRawData = new ImageStack(width, height);

            if (mode == 0) fmode = f;  // no FRC
            else if (mode == 1) fmode = 2 * (f - 1) + 1; // FRC odd frames
            else fmode = 2 * f; // FRC even frames

            imsThisRawData.addSlice(imsAllRawData.getProcessor(fmode));
            userPressedEscape = liveSRRF.calculateSRRF(imsThisRawData);

            // Check if user is cancelling calculation
            if (userPressedEscape) {
                return userPressedEscape;
            }
        }

        IJ.showStatus("Reading LiveSRRF image...");
        liveSRRF.readSRRFbuffer();

        return userPressedEscape;
    }

    private void displayImagePlus(ImageStack ims, String titleAppendix, Calibration cal, String nameLUT) {

        ImagePlus imp = new ImagePlus(imageTitle + titleAppendix, ims);
        imp.setCalibration(cal.copy());
        IJ.run(imp, "Enhance Contrast", "saturated=0.5");
        if (nameLUT.equals("ErrorMap-LUT")) applyLUT_SQUIRREL_Errors(imp);
        else if (nameLUT.equals("FRC-LUT")) applyLUT_SQUIRREL_FRC(imp);

        imp.show();
    }


}
