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

    private final int radiusCCM = 5;
    private final String LiveSRRFVersion = "v0.2";

    private float[] shiftX, shiftY;

    // Image formats
    private ImagePlus imp;
    private ImagePlus impCCM = null;

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
        CLDevice[] allDevices = liveSRRF.checkDevices();

        // Initilizaing string for device choice
        String[] deviceNames = new String[allDevices.length + 1];
        deviceNames[0] = "Default device";

        for (int i = 1; i <= allDevices.length; i++) {
            deviceNames[i] = allDevices[i - 1].getName();
        }

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

        ImageStack imsBuffer;
        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);

        ImageStack imsErrorMap = new ImageStack(width * magnification, height * magnification);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

        ImageStack imsRawData;

        int n_calculation = nframeArray.length * sensitivityArray.length * fwhmArray.length;

        IJ.log("Magnification: " + magnification);
        if (correctVibration) IJ.log("Vibration correction: on");
        else IJ.log("Vibration correction: off");

        IJ.log("Number of calculations planned: " + n_calculation);
        int r = 1;

        boolean userPressedEscape;
        ErrorMapLiveSRRF errorMapCalculator = new ErrorMapLiveSRRF();
        FloatProcessor fpErrorMap;

        int maxnFrame = nframeArray[nframeArray.length - 1];
        shiftX = new float[maxnFrame];
        shiftY = new float[maxnFrame];

        if (correctVibration) calculateShiftArray(1, maxnFrame);

        for (int thisnf : nframeArray) {
            for (int thisSensitivity : sensitivityArray) {
                for (float thisfwhm : fwhmArray) {
                    IJ.log("--------");
                    IJ.log("SRRF frame: " + r + "/" + n_calculation);
                    IJ.showProgress(r, n_calculation);

                    // Check if user is cancelling calculation
                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        liveSRRF.release();
                        return;
                    }

                    IJ.log("Number of frame for SRRF: " + thisnf);
                    IJ.log("Radius: " + thisfwhm + " pixels");
                    IJ.log("Sensitivity: " + thisSensitivity);

                    liveSRRF.initialise(width, height, magnification, thisfwhm, thisSensitivity, 1, thisnf, blockSize, null);

                    float[] shiftXtemp = new float[thisnf];
                    float[] shiftYtemp = new float[thisnf];

                    for (int i = 0; i < thisnf; i++) {
                        shiftXtemp[i] = shiftX[i];
                        shiftYtemp[i] = shiftY[i];
                    }

                    liveSRRF.loadShiftXYGPUbuffer(shiftXtemp, shiftYtemp);

                    liveSRRF.resetFramePosition();

                    for (int f = 1; f <= thisnf; f++) {
                        imp.setSlice(f);
                        imsRawData = new ImageStack(width, height);
                        imsRawData.addSlice(imp.getProcessor());
                        userPressedEscape = liveSRRF.calculateSRRF(imsRawData);

                        // Check if user is cancelling calculation
                        if (userPressedEscape) {
                            liveSRRF.release();
                            IJ.log("-------------------------------------");
                            IJ.log("Reconstruction aborted by user.");
                            return;
                        }
                    }


//                    liveSRRF.calculateSRRF(imsRawData);
                    imsBuffer = liveSRRF.readSRRFbuffer();

                    errorMapCalculator.optimise(imsBuffer.getProcessor(3), imsBuffer.getProcessor(1), magnification);
                    fpErrorMap = errorMapCalculator.calculateErrorMap();

                    if (showErrorMaps)
                        imsErrorMap.addSlice("R=" + thisfwhm + "/S=" + thisSensitivity + "/#f=" + thisnf, fpErrorMap);
                    if (calculateAVG)
                        imsSRRFavg.addSlice("R=" + thisfwhm + "/S=" + thisSensitivity + "/#f=" + thisnf, imsBuffer.getProcessor(1));
                    if (calculateSTD)
                        imsSRRFstd.addSlice("R=" + thisfwhm + "/S=" + thisSensitivity + "/#f=" + thisnf, imsBuffer.getProcessor(2));
                    imsInt.addSlice("R=" + thisfwhm + "/S=" + thisSensitivity + "/#f=" + thisnf, imsBuffer.getProcessor(3));

                    r++;
                }
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

        IJ.log("-------------------------------------");
        IJ.log("RAM used: " + IJ.freeMemory());
        IJ.log("Bye-bye !");

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
    private void calculateShiftArray(int indexStart, int nFrameForSRRF) {

        imp.setSlice(indexStart);
        ImageProcessor ipRef = imp.getProcessor().duplicate();
        ImageProcessor ipData;

        for (int s = 0; s < nFrameForSRRF; s++) {

            // Check if user is cancelling calculation
            IJ.showProgress(s, nFrameForSRRF);
            if (IJ.escapePressed()) {
                liveSRRF.release();
                IJ.resetEscape();
                IJ.log("-------------------------------------");
                IJ.log("Reconstruction aborted by user.");
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


}
