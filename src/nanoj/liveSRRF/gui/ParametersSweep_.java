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
            showErrorMaps,
            showRSC,
            calculateFRC,
            calculateRSE,
            calculateRSP,
            doErrorMapping,
            fixSigma,
            correctVibration,
            cropBorder;

    private float[] fwhmArray;
    private int[] sensitivityArray,
            nframeArray;

    private float fixedSigma;

    private final String LiveSRRFVersion = "v1.1";
    private float[] shiftX, shiftY;

    private String chosenTemporalAnalysis;

    // Advanced formats
    private final NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private final NanoJProfiler prof = new NanoJProfiler();
//    private liveSRRF_CL liveSRRF;

    public void run(String arg) {


        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
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

        // Initialise the liveSRRF engine (creates a context really)
        liveSRRF_CL liveSRRF = new liveSRRF_CL();

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        String[] temporalAnalysis = {"AVG","STD","Both AVG and STD"};


        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("liveSRRF - Parameters sweep " + LiveSRRFVersion);
        gd.addMessage("-=-= liveSRRF reconstruction =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addChoice("Temporal analysis", temporalAnalysis, prefs.get("chosenTemporalAnalysis", temporalAnalysis[2]));
        gd.addCheckbox("Correct vibration (default: off)", prefs.get("correctVibration", false));
        gd.addCheckbox("Show all reconstructions (default: on)", prefs.get("showRecons", true));

        gd.addMessage("-=-= Sweeping liveSRRF parameters =-=-\n", headerFont);
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
            IJ.log("Error mapping paramaters:");
            if (fixSigma) IJ.log("Sigma is fixed to " + fixedSigma + " pixels.");
            else IJ.log("Sigma is optimised for each reconstructions.");
            if (cropBorder) IJ.log("The borders will be cropped by "+cropSize+" pixels");

        }


        int n_calculation = nframeArray.length * sensitivityArray.length * fwhmArray.length;
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
        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);


        // Initialising the Error maps variables
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF(imp, magnification, fixSigma, fixedSigma, cropBorder, cropSize);

        float[] pixelsRMSEavg = null;
        float[] pixelsPPMCCavg = null;
        float[] pixelsRMSEstd = null;
        float[] pixelsPPMCCstd = null;

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
        ImageStack imsErrorMapSTD = new ImageStack(widthSquirrel, heightSquirrel);

        ImageStack imsRSCavg = new ImageStack(widthSquirrel, heightSquirrel);
        ImageStack imsRSCstd = new ImageStack(widthSquirrel, heightSquirrel);



        // These images are small so it doesn't matter if they are initialised in any case
        ImageStack imsRMSEavg = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCCavg = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsRMSEstd = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsPPMCCstd = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);


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
        FloatProcessor fpOddSTD;
        FloatProcessor fpEvenSTD;

        float[] pixelsFRCresolutionAVG = null;
        float[] pixelsFRCresolutionSTD = null;

        // These images are small so it doesn't matter if they are initialised in any case
        ImageStack imsFRCresolutionAVG = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);
        ImageStack imsFRCresolutionSTD = new ImageStack(fwhmArray.length, sensitivityArray.length, nframeArray.length);

        // Set the number of reconstuctions
        int r = 0;
        boolean userPressedEscape;

        for (int nfi = 0; nfi < nframeArray.length; nfi++) {

            if (doErrorMapping) {
                pixelsRMSEavg = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsPPMCCavg = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsRMSEstd = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsPPMCCstd = new float[(fwhmArray.length) * (sensitivityArray.length)];
            }

            if (calculateFRC) {
                pixelsFRCresolutionAVG = new float[(fwhmArray.length) * (sensitivityArray.length)];
                pixelsFRCresolutionSTD = new float[(fwhmArray.length) * (sensitivityArray.length)];
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
//                IJ.log("ShiftX="+shiftXtemp[i]);
                    shiftYtempOdd[i] = shiftY[2 * i + 1];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }

                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXtempEven[i] = shiftX[2 * i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
                    shiftYtempEven[i] = shiftY[2 * i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }

            } else {
                for (int i = 0; i < nframeArray[nfi]; i++) {
                    shiftXtemp[i] = shiftX[i];
//                IJ.log("ShiftX="+shiftXtemp[i]);
                    shiftYtemp[i] = shiftY[i];
//                IJ.log("ShiftY="+shiftYtemp[i]);
                }
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

                        // Calculate and get the reconstruction from the odd frames // TODO: add options to do intensity weighting or MPcorrection
                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, true, null, 0,1);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXtempOdd, shiftYtempOdd);

                        userPressedEscape = calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 1, liveSRRF);
                        if (userPressedEscape) {
                            liveSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }

                        imsBuffer = liveSRRF.imsSRRF;

                        fpOddAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                        fpOddSTD = imsBuffer.getProcessor(2).convertToFloatProcessor();

                        // Calculate and get the reconstruction from the even frames // TODO: add options to do intensity weighting or MPcorrection
                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, true, null, 0, 1);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXtempEven, shiftYtempEven);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 2, liveSRRF);
                        imsBuffer = liveSRRF.imsSRRF;

                        if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                            fpEvenAVG = imsBuffer.getProcessor(1).convertToFloatProcessor();
                            pixelsFRCresolutionAVG[fwhmArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddAVG, fpEvenAVG, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        }
                        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                            fpEvenSTD = imsBuffer.getProcessor(2).convertToFloatProcessor();
                            pixelsFRCresolutionSTD[fwhmArray.length * si + fi] = (float) cal.pixelHeight * (float) frcCalculator.calculateFireNumber(fpOddSTD, fpEvenSTD, FRC.ThresholdMethod.FIXED_1_OVER_7);
                        }


                    } else {
                        // TODO: add options to do intensity weighting or MPcorrection
                        liveSRRF.initialise(width, height, magnification, fwhmArray[fi], sensitivityArray[si], 1, nframeArray[nfi], blockSize, null, true, true, null, 0, 1);
                        liveSRRF.resetFramePosition();
                        liveSRRF.loadDriftXYGPUbuffer(shiftXtemp, shiftYtemp);

                        calculateLiveSRRFsingleframeLoad(imsAllRawData, nfi, 0, liveSRRF);
                        imsBuffer = liveSRRF.imsSRRF;
                    }


                    // Collate the reconstructions
                    if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2]))
                        imsSRRFavg.addSlice(label, imsBuffer.getProcessor(1));
                    if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2]))
                        imsSRRFstd.addSlice(label, imsBuffer.getProcessor(2));

                    // Interpolated image
                    imsInt.addSlice(label, imsBuffer.getProcessor(3));

                    // Error mapping
                    if (doErrorMapping) {

                        if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                            // Error maps for AVG
                            IJ.showStatus("Optimising Sigma...");
                            errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(1));

                            IJ.showStatus("Calculating error map...");
                            errorMapCalculator.calculateErrorMap();
                            pixelsRMSEavg[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                            pixelsPPMCCavg[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                            if (showErrorMaps) imsErrorMapAVG.addSlice(label, errorMapCalculator.fpErrorMap);
                            if (showRSC) imsRSCavg.addSlice(label, errorMapCalculator.fpSRC);

                        }

                        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                            // Error maps for STD
                            IJ.showStatus("Optimising Sigma...");
                            errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(2));

                            IJ.showStatus("Calculating error map...");
                            errorMapCalculator.calculateErrorMap();
                            pixelsRMSEstd[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalRMSE;
                            pixelsPPMCCstd[fwhmArray.length * si + fi] = (float) errorMapCalculator.globalPPMCC;

                            if (showErrorMaps) imsErrorMapSTD.addSlice(label, errorMapCalculator.fpErrorMap);
                            if (showRSC) imsRSCstd.addSlice(label, errorMapCalculator.fpSRC);
                        }
                    }


                    // Increment the reconstruction counter
                    r++;
                }
            }

            if (calculateRSE) {
                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsRMSEavg.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSEavg), nfi + 1);
                    imsRMSEavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }

                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsRMSEstd.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsRMSEstd), nfi + 1);
                    imsRMSEstd.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }

            if (calculateRSP) {
                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsPPMCCavg.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCCavg), nfi + 1);
                    imsPPMCCavg.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsPPMCCstd.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsPPMCCstd), nfi + 1);
                    imsPPMCCstd.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }

            if (calculateFRC) {
                if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsFRCresolutionAVG.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionAVG), nfi + 1);
                    imsFRCresolutionAVG.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
                if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                    imsFRCresolutionSTD.setProcessor(new FloatProcessor(fwhmArray.length, sensitivityArray.length, pixelsFRCresolutionSTD), nfi + 1);
                    imsFRCresolutionSTD.setSliceLabel("#fr=" + nframeArray[nfi], nfi + 1);
                }
            }

        }

        // Release the GPU
        liveSRRF.release();

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


        // liveSRRF (AVG) reconstruction
        if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG)", imsSRRFavg);
            impSRRFavg.setCalibration(cal);
            IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
            impSRRFavg.show();
        }

        // liveSRRF (STD) reconstruction
        if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD)", imsSRRFstd);
            impSRRFstd.setCalibration(cal);
            IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
            impSRRFstd.show();
        }

        // Interpolated image
        ImagePlus impInt = new ImagePlus(imp.getTitle() + " - Interpolated image", imsInt);
        impInt.setCalibration(cal);
        IJ.run(impInt, "Enhance Contrast", "saturated=0.5");
        impInt.show();

        // Showing error maps
        if (showErrorMaps) {
            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impErrorMapAVG = new ImagePlus(imp.getTitle() + " - ErrorMap (AVG)", imsErrorMapAVG);
                impErrorMapAVG.setCalibration(cal);
                IJ.run(impErrorMapAVG, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impErrorMapAVG);
                impErrorMapAVG.show();
            }
            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impErrorMapSTD = new ImagePlus(imp.getTitle() + " - ErrorMap (STD)", imsErrorMapSTD);
                impErrorMapSTD.setCalibration(cal);
                IJ.run(impErrorMapSTD, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impErrorMapSTD);
                impErrorMapSTD.show();
            }
        }

        if (showRSC) {
            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impRSCavg = new ImagePlus(imp.getTitle() + " - rescaled liveSRRF (AVG)", imsRSCavg);
                impRSCavg.setCalibration(cal);
                IJ.run(impRSCavg, "Enhance Contrast", "saturated=0.5");
                impRSCavg.show();
            }
            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impRSCstd = new ImagePlus(imp.getTitle() + " - rescaled liveSRRF (STD)", imsRSCstd);
                impRSCstd.setCalibration(cal);
                IJ.run(impRSCstd, "Enhance Contrast", "saturated=0.5");
                impRSCstd.show();
            }
        }

        if (calculateRSE) {
            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impRMSEavg = new ImagePlus(imp.getTitle() + " - RMSE sweep map (AVG)", imsRMSEavg);
                impRMSEavg.setCalibration(sweepMapCalib.copy());
                IJ.run(impRMSEavg, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impRMSEavg);
                impRMSEavg.show();
            }
            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impRMSEstd = new ImagePlus(imp.getTitle() + " - RMSE sweep map (STD)", imsRMSEstd);
                impRMSEstd.setCalibration(sweepMapCalib.copy());
                IJ.run(impRMSEstd, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impRMSEstd);
                impRMSEstd.show();
            }
        }

        if (calculateRSP) {
            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impPPMCCavg = new ImagePlus(imp.getTitle() + " - RSP sweep map (AVG)", imsPPMCCavg);
                impPPMCCavg.setCalibration(sweepMapCalib.copy());
                IJ.run(impPPMCCavg, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impPPMCCavg);
                impPPMCCavg.show();
            }
            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impPPMCCstd = new ImagePlus(imp.getTitle() + " - RSP sweep map (STD)", imsPPMCCstd);
                impPPMCCstd.setCalibration(sweepMapCalib.copy());
                IJ.run(impPPMCCstd, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_Errors(impPPMCCstd);
                impPPMCCstd.show();
            }
        }

        if (calculateFRC) {
            if (chosenTemporalAnalysis.equals(temporalAnalysis[0]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impFRCresolutionAVG = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (AVG)", imsFRCresolutionAVG);
                impFRCresolutionAVG.setCalibration(sweepMapCalib.copy());
                IJ.run(impFRCresolutionAVG, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_FRC(impFRCresolutionAVG);
                impFRCresolutionAVG.show();
            }
            if (chosenTemporalAnalysis.equals(temporalAnalysis[1]) || chosenTemporalAnalysis.equals(temporalAnalysis[2])) {
                ImagePlus impFRCresolutionSTD = new ImagePlus(imp.getTitle() + " - FRC resolution sweep map (STD)", imsFRCresolutionSTD);
                impFRCresolutionSTD.setCalibration(sweepMapCalib.copy());
                IJ.run(impFRCresolutionSTD, "Enhance Contrast", "saturated=0.5");
                applyLUT_SQUIRREL_FRC(impFRCresolutionSTD);
                impFRCresolutionSTD.show();
            }
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
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    //    --- Grab settings ---
    private boolean grabSettings(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();
        chosenTemporalAnalysis = gd.getNextChoice();

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
        // TODO: ABORT IF NOT enough frames in original stack for even 1 analysis

        if (n_nfToUse < 0){
            return false;
        }

        nframeArray = new int[n_nfToUse];
        for (int i = 0; i < n_nfToUse; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }


        prefs.set("magnification", magnification);
        prefs.set("chosenTemporalAnalysis", chosenTemporalAnalysis);
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

        prefs.set("blockSize", blockSize);


        prefs.save();

        return true;
    }


    private boolean calculateLiveSRRFsingleframeLoad(ImageStack imsAllRawData, int nf, int mode, liveSRRF_CL liveSRRF) {

        ImageStack imsThisRawData;
        boolean userPressedEscape = false;
        int fmode;

        IJ.showStatus("Calculating liveSRRF image...");
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

        IJ.showStatus("Reading liveSRRF image...");
        liveSRRF.readSRRFbuffer();

        return userPressedEscape;

    }


}
