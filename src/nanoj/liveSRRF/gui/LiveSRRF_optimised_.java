package nanoj.liveSRRF.gui;

import com.jogamp.opencl.CLDevice;
import ij.*;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import nanoj.core.java.io.zip.SaveFileInZip;
import nanoj.core.java.io.zip.virtualStacks.FullFramesVirtualStack;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.core2.NanoJUsageTracker;
import nanoj.liveSRRF.XYShiftCalculator;
import org.python.modules.math;
import nanoj.liveSRRF.LiveSRRF_CL;

import java.awt.*;
import java.io.IOException;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

import static java.lang.Math.min;

public class LiveSRRF_optimised_ implements PlugIn {

    // Basic formats
    private int magnification,
            nFrameForSRRFtoUse,
            nFrameForSRRFfromUser,
            sensitivity,
            frameGapToUse,
            frameGapFromUser,
            nFrameOnGPU,
            nSRRFframe,
            nSlices,
            width, height,
            widthS, heightS,
            blockSize,
            nGPUloadPerSRRFframe,
            gradMag,
            nRecons,
            nPlanes,
            nPlanesM;

    private float fwhm, maxMemoryGPU, axialOffset;

    private final float maxMemoryRAMij = (float) IJ.maxMemory()/1e6f; // maximum RAM set for Fiji in MB

    private boolean correctVibration,
            calculateAVG,
            calculateVAR,
            calculateTAC2,
            getInterpolatedImage,
            writeToDiskToUse,
            writeToDiskTemp,
            doRollingAnalysis,
            previousWriteToDisk = false,
            previousAdvSettings = false,
            writeSuggestOKeyed = false,
            doMPmapCorrection,
            showImStabPlot,
            intWeighting,
            showGradient,
            showIntGradient,
            do3DSRRF,
            DEBUG = false;

    private final boolean enable3DSRRF = true;
    private final String eSRRFVersion = "v1.4.0";

    private String pathToDisk = "",
            fileName,
            chosenDeviceName;

    private String[] deviceNames;

    private float[] shiftX, shiftY;

    // Image formats
    private ImagePlus imp,
            impSRRFavg,
            impSRRFvar,
            impSRRFtac2,
            impRawInterpolated;


    // Advanced formats
    private final NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private final NanoJProfiler prof = new NanoJProfiler();
    private LiveSRRF_CL eSRRF;

    private SaveFileInZip saveFileInZip;
    private CLDevice chosenDevice = null;

    // Tracker
    private final String user = "FijiUpdater";
    private final String version = "20180809-" + user;
    private final NanoJUsageTracker tracker = new NanoJUsageTracker("NanoJ-LiveSRRF", version, "UA-61590656-4");


    @Override
    public void run(String arg) {

        tracker.logUsage("LiveSRRF_");

        // Get raw data
        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double pixelSize = imageInfo.pixelHeight; //um if it is
        String pixelSizeUnit = imageInfo.unit;

        nSlices = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();

        // initialization of advanced settings (in case advanced settings are not selected)
        chosenDeviceName = prefs.get("chosenDeviceName", "Default device");
        maxMemoryGPU = prefs.get("maxMemoryGPU", 500);
        blockSize = (int) prefs.get("blockSize", 20000);
        writeToDiskToUse = false;
        intWeighting = prefs.get("intWeighting", true);
        doMPmapCorrection = prefs.get("doMPmapCorrection", true);

//        // Close the log: //TODO: this seems to sometimes make a different whether OpenCL runs or not (observation from Nvidia 1050, not repeatable)
        // TODO: but this causes the GUI to pop up when using liveSRRF in a macro
//        if (IJ.getLog() != null) {
//            selectWindow("Log");
//            IJ.run("Close");
//        }

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("eSRRF " + eSRRFVersion); //TODO: what's taking so long on Wolverine?
        IJ.log("Max RAM available: "+ (float) Math.round(maxMemoryRAMij*100)/100 + " MB");

        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));
        IJ.log("Pixel size: "+pixelSize+" "+pixelSizeUnit);

        // Initialize the liveSRRF class and check the devices
        eSRRF = new LiveSRRF_CL(DEBUG);
        eSRRF.checkDevices();
        CLDevice[] allDevices = LiveSRRF_CL.allCLdevices;
        gradMag = eSRRF.gradientMag;
        nRecons = eSRRF.nReconstructions;

        // Initializing string for device choice
        deviceNames = new String[allDevices.length + 1];
        deviceNames[0] = "Default device";

        for (int i = 1; i <= allDevices.length; i++) {
            deviceNames[i] = allDevices[i - 1].getName();
        }

        // Main GUI ----
        boolean cancelled = mainGUI();
        if (cancelled) {
//            eSRRF.release(); // no need te release contex prior to initilisation of liveSRRF
            IJ.log("Cancelled by user.");
            return;
        }

        // Check if something has gone wrong with the memory
        if (!calculatenFrameOnGPU()){
//            eSRRF.release(); // no need te release contex prior to initilisation of liveSRRF
            IJ.log("Problems with memory. Check advanced settings or else...");
            return;
        }

        // Get chosen device
        if (chosenDeviceName == null) chosenDeviceName = "Default device";
        for (CLDevice thisDevice : allDevices) {
            if (chosenDeviceName.equals(thisDevice.getName())) chosenDevice = thisDevice;
        }

        // ---- Getting calibration data from the NanoJ table ----
        String calibTablePath3DSRRF = null;
        if (do3DSRRF){
            IJ.log("Getting 3D-SRRF calibration file...");
            calibTablePath3DSRRF = IJ.getFilePath("Choose 3D-SRRF calibration table to load...");
            if (calibTablePath3DSRRF == null) return;
        }

        // Save last user entries
        savePreferences();

        IJ.log("-------------------------------------");
        IJ.log("Parameters:");
        IJ.log("Magnification: " + magnification);
        IJ.log("Radius: " + fwhm + " pixels");
        IJ.log("Sensitivity: " + sensitivity);
        IJ.log("# frames for SRRF: " + nFrameForSRRFtoUse);
        IJ.log("# frames gap: " + frameGapToUse);
        IJ.log("# frames on device: " + nFrameOnGPU);
        IJ.log("Estimated device memory usage: " + (float) Math.ceil(predictMemoryUsed(nFrameOnGPU)[0] * 100)/100 + " MB");
        IJ.log("Estimated RAM usage: " + (float) Math.ceil(predictMemoryUsed(nFrameOnGPU)[1] * 100)/100 + " MB");
        IJ.log("Running on: " + chosenDeviceName);
        IJ.log("# reconstructed frames: " + nSRRFframe);
        IJ.log("# device load / SRRF frame: " + nGPUloadPerSRRFframe);
        IJ.log("Macro-pixel artefact removal: " + doMPmapCorrection);
        IJ.log("Vibration correction: "+ correctVibration);

        // Initialize ZipSaver in case we're writing to disk
        if (writeToDiskToUse) {
            fileName = pathToDisk + imp.getTitle() + " - eSRRF.zip";
            IJ.log("Write to disk enabled");
            IJ.log(fileName);

            try {
                saveFileInZip = new SaveFileInZip(fileName, false);
            } catch (IOException e) {
                IJ.error("Whoops, it seems that " + fileName + " doesn't exist. At least there's ice-cream...");
                e.printStackTrace();
            }
        }

        // Initialize variables
        int indexStartSRRFframe;
        int nFrameToLoad;

        ImageStack imsAllRawData = imp.getImageStack();

        eSRRF.initialise(width, height, magnification, fwhm, sensitivity, nFrameOnGPU, nFrameForSRRFtoUse, blockSize, chosenDevice, intWeighting, doMPmapCorrection, calibTablePath3DSRRF, axialOffset, (float) pixelSize);
        widthS = eSRRF.widthS;
        heightS = eSRRF.heightS;
        nPlanes = eSRRF.nPlanes;
        nPlanesM = eSRRF.nPlanesM;
//        IJ.log("WidthS/HeightS: "+widthS+"/"+heightS);
        IJ.log("Number of planes: "+nPlanes);

        shiftX = new float[nFrameForSRRFtoUse];
        shiftY = new float[nFrameForSRRFtoUse];

        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp, prof);
        ImageStack imsBuffer;

        ImageStack imsRawData;
        ImageStack imsSRRFavg = new ImageStack(widthS * magnification, heightS * magnification);
        ImageStack imsSRRFvar = new ImageStack(widthS * magnification, heightS * magnification);
        ImageStack imsSRRFtac2 = new ImageStack(widthS * magnification, heightS * magnification);
        ImageStack imsRawInterpolated = new ImageStack(widthS * magnification, heightS * magnification);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

        if (DEBUG) IJ.log("Distance units: "+cal.getUnit());

        boolean userPressedEscape;
        long loopStart = System.nanoTime();

        // For measuring the remaining time
        double singleLoadTime;
        double remainingTime;
        int currentNLoad;

        // Start looping trough SRRF frames --------------------------------------------
        for (int r = 1; r <= nSRRFframe; r++) {
            eSRRF.resetFramePosition();  // resets only the global SRRF counter

            IJ.log("--------");
            IJ.log("SRRF frame: " + r + "/" + nSRRFframe);
            indexStartSRRFframe = (r - 1) * frameGapToUse + 1;
            IJ.log("Stack index start: " + indexStartSRRFframe);

            if (correctVibration) {
                shiftCalculator.calculateShiftArray(indexStartSRRFframe, nFrameForSRRFtoUse);
                shiftX = shiftCalculator.shiftX;
                shiftY = shiftCalculator.shiftY;

                if (showImStabPlot) {
                    // Scatter plot
                    Plot scatterPlot = new Plot("x/y scatter plot (spectrum LUT-coded) #" + r, "x (pixels)", "y (pixels)");
                    for (int i = 0; i < nFrameForSRRFtoUse; i++) {
                        float[] x_temp = new float[1];
                        float[] y_temp = new float[1];
                        x_temp[0] = shiftX[i];
                        y_temp[0] = shiftY[i];
                        scatterPlot.setColor(Color.getHSBColor(i / (float) nSlices, 1f, 1f)); // this corresponds to the spectrum LUT
                        scatterPlot.addPoints(x_temp, y_temp, Plot.CROSS);
                    }
                    scatterPlot.show();
                    scatterPlot.setLimitsToFit(true);
                }
            }

            eSRRF.loadDriftXYGPUbuffer(shiftX, shiftY);

            IJ.showProgress(r - 1, nSRRFframe);

            for (int l = 0; l < nGPUloadPerSRRFframe; l++) {

                nFrameToLoad = min(nFrameOnGPU, nFrameForSRRFtoUse - nFrameOnGPU * l);
                if (DEBUG) {
                    IJ.log("----------------------------");
                    IJ.log("Number of frames to load: " + nFrameToLoad);
                    IJ.log("Index start: " + (indexStartSRRFframe + nFrameOnGPU * l));
                }

                imsRawData = new ImageStack(width, height);
                for (int f = 0; f < nFrameToLoad; f++) {
                    imsRawData.addSlice(imsAllRawData.getProcessor(indexStartSRRFframe + nFrameOnGPU * l + f));
                }

                userPressedEscape = eSRRF.calculateSRRF(imsRawData); // resets the local GPU load frame counter

                // Check if user is cancelling calculation
                if (userPressedEscape) {
                    eSRRF.release();
                    IJ.log("-------------------------------------");
                    IJ.log("Reconstruction aborted by user.");
                    return;
                }

                // Estimating the remaining time
                currentNLoad = nGPUloadPerSRRFframe*(r-1) + l+1;
                singleLoadTime = ((System.nanoTime() - loopStart) / (float) currentNLoad) / 1e9;
                remainingTime = singleLoadTime * (nSRRFframe * nGPUloadPerSRRFframe - currentNLoad);
                IJ.showStatus("eSRRF - Remaining time: " + timeToString(remainingTime));
            }

            eSRRF.readSRRFbuffer();
            imsBuffer = eSRRF.imsSRRF;

            if (writeToDiskToUse) {//  TODO: this will not work in 3D
                try {
                    if (calculateAVG) {
                        impTemp = new ImagePlus("\"eSRRF (AVG) - frame=" + r + ".tif", imsBuffer.getProcessor(1));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("eSRRF (AVG) - frame=" + r, impTemp);
                    }
                    if (calculateVAR) {
                        impTemp = new ImagePlus("\"eSRRF (VAR) - frame=" + r + ".tif", imsBuffer.getProcessor(2));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("eSRRF (VAR) - frame=" + r, impTemp);
                    }
                    if (calculateTAC2) {
                        impTemp = new ImagePlus("\"eSRRF (TAC2) - frame=" + r + ".tif", imsBuffer.getProcessor(3));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("eSRRF (TAC2) - frame=" + r, impTemp);
                    }
                    if (getInterpolatedImage) {
                        impTemp = new ImagePlus("\"eSRRF (INT) - frame=" + r + ".tif", imsBuffer.getProcessor(4));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("eSRRF (INT) - frame=" + r, impTemp);
                    }
                } catch (IOException e) {
                    IJ.error("Whoops, it seems that there was a problem with saving to disk... (insert sad face here).");
                    e.printStackTrace();
                }

            } else {

                    if (calculateAVG) {
                        for (int p = 0; p < nPlanesM; p++) imsSRRFavg.addSlice(imsBuffer.getProcessor(p+1));
                    }
                    if (calculateVAR) {
                        for (int p = 0; p < nPlanesM; p++) imsSRRFvar.addSlice(imsBuffer.getProcessor(p + nPlanesM+1));
                    }
                if (calculateTAC2) {
                    for (int p = 0; p < nPlanesM; p++) imsSRRFtac2.addSlice(imsBuffer.getProcessor(p + 2*nPlanesM+1));
                }
                    if (getInterpolatedImage) {
                        for (int p = 0; p < nPlanesM; p++) imsRawInterpolated.addSlice(imsBuffer.getProcessor(p + 3*nPlanesM+1));
                    }
                }

            IJ.log("RAM used: " + IJ.freeMemory());
        }
        // End looping trough SRRF frames --------------------------------------------

        if (showGradient) { // TODO: small bug leading to not remembering choices of showGradient
            // Read off the gradient (single frame, which one needs to be checked)
            ImageStack imsGradient = eSRRF.readGradientBuffers(false);
            ImagePlus impGradients = new ImagePlus("Gradients", imsGradient);
            impGradients.show();
        }

        if (showIntGradient) {
            ImageStack imsGradientInt = eSRRF.readGradientBuffers(true);
            ImagePlus impGradientsInt = new ImagePlus("Interpolated gradients", imsGradientInt);
            impGradientsInt.show();
        }

        // ----- DEBUG -----
        if (DEBUG) {
            if (doMPmapCorrection) {
                ImageStack imsMPmap = eSRRF.readMPmaps();
                ImagePlus impMPmap = new ImagePlus("MP map", imsMPmap);
                impMPmap.show();
            }
            if (do3DSRRF) {
                ImageStack imsAlignedPixels = eSRRF.readAlignedPixels();
                ImagePlus impAlignedPixels = new ImagePlus("Aligned pixels", imsAlignedPixels);
                impAlignedPixels.show();
            }
        }

        // Release the GPU
        eSRRF.release();

        IJ.log("-------------------------------------");
        // Close the ZipSaver & display stack as virtual stacks
        if (writeToDiskToUse) {
            IJ.log("Results displayed as virtual stacks.");

            try {
                saveFileInZip.close();
            } catch (IOException e) {
                IJ.error("Error closing the ZipSaver !!!!!");
                e.printStackTrace();
            }

            try {
                if (calculateAVG) {
                    FullFramesVirtualStack vsimsSRRFavg = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFavg.getSize(); f > 0; f--) {
                        if (!vsimsSRRFavg.getSliceLabel(f).contains("AVG")) vsimsSRRFavg.deleteSlice(f);
                    }

                    impSRRFavg = new ImagePlus(imp.getTitle() + " - eSRRF (AVG)", vsimsSRRFavg);
                    impSRRFavg.setCalibration(cal);
                    IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                    impSRRFavg.show();
                }

                if (calculateVAR) {
                    FullFramesVirtualStack vsimsSRRFvar = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFvar.getSize(); f > 0; f--) {
                        if (!vsimsSRRFvar.getSliceLabel(f).contains("VAR")) vsimsSRRFvar.deleteSlice(f);
                    }

                    impSRRFvar = new ImagePlus(imp.getTitle() + " - eSRRF (VAR)", vsimsSRRFvar);
                    impSRRFvar.setCalibration(cal);
                    IJ.run(impSRRFvar, "Enhance Contrast", "saturated=0.5");
                    impSRRFvar.show();
                }

                if (calculateTAC2) {
                    FullFramesVirtualStack vsimsSRRFtac2 = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFtac2.getSize(); f > 0; f--) {
                        if (!vsimsSRRFtac2.getSliceLabel(f).contains("TAC2")) vsimsSRRFtac2.deleteSlice(f);
                    }

                    impSRRFtac2 = new ImagePlus(imp.getTitle() + " - eSRRF (TAC2)", vsimsSRRFtac2);
                    impSRRFtac2.setCalibration(cal);
                    IJ.run(impSRRFtac2, "Enhance Contrast", "saturated=0.5");
                    impSRRFtac2.show();
                }

                if (getInterpolatedImage) {
                    FullFramesVirtualStack vsimsRawInterpolated = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsRawInterpolated.getSize(); f > 0; f--) {
                        if (!vsimsRawInterpolated.getSliceLabel(f).contains("INT")) vsimsRawInterpolated.deleteSlice(f);
                    }

                    impRawInterpolated = new ImagePlus(imp.getTitle() + " - interpolated image", vsimsRawInterpolated);
                    impRawInterpolated.setCalibration(cal);
                    IJ.run(impRawInterpolated, "Enhance Contrast", "saturated=0.5");
                    impRawInterpolated.show();
                }

            } catch (IOException e) {
                IJ.error("Whoops, it seems that there was a problem with opening from disk... (insert sad face here).");
                e.printStackTrace();
            }


        } else {

            //Display results
            if (calculateAVG) {
                impSRRFavg = new ImagePlus(imp.getTitle() + " - eSRRF (AVG)", imsSRRFavg);
                impSRRFavg.setCalibration(cal);
                IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                impSRRFavg.show(); //TODO: add z-step size for 3D
            }

            if (calculateVAR) {
                impSRRFvar = new ImagePlus(imp.getTitle() + " - eSRRF (VAR)", imsSRRFvar);
                impSRRFvar.setCalibration(cal);
                IJ.run(impSRRFvar, "Enhance Contrast", "saturated=0.5");
                impSRRFvar.show();
            }

            if (calculateTAC2) {
                impSRRFtac2 = new ImagePlus(imp.getTitle() + " - eSRRF (TAC2)", imsSRRFtac2);
                impSRRFtac2.setCalibration(cal);
                IJ.run(impSRRFtac2, "Enhance Contrast", "saturated=0.5");
                impSRRFtac2.show();
            }

            if (getInterpolatedImage) {
                impRawInterpolated = new ImagePlus(imp.getTitle() + " - interpolated image", imsRawInterpolated);
                impRawInterpolated.setCalibration(cal);
                IJ.run(impRawInterpolated, "Enhance Contrast", "saturated=0.5");
                impRawInterpolated.show();
            }
        }


        // Bye-bye and report
        IJ.log("Memory usage: " + IJ.freeMemory());  // this also runs the garbage collector
        IJ.log("Thank you for your custom on this beautiful day !");
        now = LocalDateTime.now();
        IJ.log(now.format(formatter));

        long executionEndTime = System.nanoTime();
        double executionTime = (executionEndTime - loopStart)/1e9;
        IJ.log("Execution time: " + timeToString(executionTime));

        // Run garbage collector
        System.gc();
        IJ.run("Cascade", "");

        // ----- DEBUG -----
        if (DEBUG) {
            IJ.log("-------------------------------------");
            IJ.log(prof.report());
        }
    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    // -- Main GUI --
    private boolean mainGUI() {
        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("eSRRF " + eSRRFVersion);
        gd.addMessage("-=-= SRRF parameters =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("Radius (pixels, default: 1.5)", prefs.get("fwhm", (float) 1.5), 2);
        gd.addNumericField("Sensitivity (default: 1)", prefs.get("sensitivity", 1), 0);
        gd.addNumericField("# frames for SRRF (0 = auto)", prefs.get("nFrameForSRRFfromUser", 0), 0);
        gd.addCheckbox("Vibration correction", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("VAR reconstruction (default: off)", prefs.get("calculateVAR", false));
        gd.addCheckbox("TAC2 reconstruction (default: off)", prefs.get("calculateTAC2", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("getInterpolatedImage", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addCheckbox("Perform rolling analysis (default: off)", prefs.get("doRollingAnalysis", false));
        gd.addNumericField("# frame gap between SR frame (0 = auto)", prefs.get("frameGapFromUser", 0), 0);

        gd.addMessage("-=-= Advanced settings =-=-\n", headerFont);
        gd.addCheckbox("Show advanced settings", false);
        if (enable3DSRRF) {
            gd.addCheckbox("3D-SRRF from MFM data", prefs.get("do3DSRRF", false));
            gd.addNumericField("Axial offset (in nm): ",400, 2);
        }

        gd.addHelp("https://www.youtube.com/watch?v=PJQVlVHsFF8"); // If you're hooked on a feeling

        MyDialogListenerMainGUI dl = new MyDialogListenerMainGUI(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        grabSettingsMainGUI(gd); // important to make sure that settings are grabbed properly even when no parameters is changed

        gd.showDialog();

        // If the GUI was cancelled
        return gd.wasCanceled();
    }

    // --- Main GUI Dialog listener ---
    class MyDialogListenerMainGUI implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {
            return grabSettingsMainGUI(gd);
        }
    }

    //    --- Grab settings from main GUI ---
    private boolean grabSettingsMainGUI(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();
        fwhm = (float) gd.getNextNumber();
        sensitivity = (int) gd.getNextNumber();
        nFrameForSRRFfromUser = (int) gd.getNextNumber();

        correctVibration = gd.getNextBoolean();

        calculateAVG = gd.getNextBoolean();
        calculateVAR = gd.getNextBoolean();
        calculateTAC2 = gd.getNextBoolean();
        getInterpolatedImage = gd.getNextBoolean();

        doRollingAnalysis = gd.getNextBoolean();
        frameGapFromUser = (int) gd.getNextNumber();

        // Calculate the parameters based on user inputs
        if (nFrameForSRRFfromUser == 0) nFrameForSRRFtoUse = nSlices;
        else nFrameForSRRFtoUse = min(nSlices, nFrameForSRRFfromUser);

        if (frameGapFromUser == 0 || !doRollingAnalysis) frameGapToUse = nFrameForSRRFtoUse;
        else frameGapToUse = frameGapFromUser;

        nSRRFframe = (int) ((float) (nSlices - nFrameForSRRFtoUse) / (float) frameGapToUse) + 1;

        boolean goodToGo = calculatenFrameOnGPU();

        boolean showAdvancedSettings = gd.getNextBoolean();
        if (showAdvancedSettings && !previousAdvSettings) advancedSettingsGUI();
        previousAdvSettings = showAdvancedSettings;

        if (enable3DSRRF) {
            do3DSRRF = gd.getNextBoolean();
            axialOffset = (float) gd.getNextNumber();
        }
        else {
            do3DSRRF = false;
            axialOffset = 0; // this doesn't matter if 3D is not performed
        }

        nGPUloadPerSRRFframe = (int) math.ceil((float) nFrameForSRRFtoUse / (float) nFrameOnGPU);

        return goodToGo;
    }

    // -- Advanced settings GUI --
    private void advancedSettingsGUI() {

        double[] memUsed = predictMemoryUsed(1);

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        GenericDialog gd = new GenericDialog("eSRRF - Advanced settings");

        gd.addMessage("-=-= GPU/CPU processing =-=-\n", headerFont);
        gd.addChoice("Processing device", deviceNames, prefs.get("chosenDeviceName", "Default device"));
        gd.addNumericField("Maximum amount of memory on device (MB, default: 1000)", prefs.get("maxMemoryGPU", 500), 2);
        gd.addMessage("Minimum device memory necessary: " + (float) Math.ceil(memUsed[0] * 100) / 100 + "MB\n");

        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use\n" +
                "more resources and may slow down your computer.");

        gd.addMessage("-=-= Write to disk =-=-\n", headerFont);
        gd.addCheckbox("Directly write to disk (default: off)", writeToDiskToUse);
        gd.addMessage("Writing directly to disk will slow down the reconstruction but \n" +
                "will allow for long time courses to be reconstructed without\n" +
                "exceeding RAM capacity.");

        gd.addMessage("-=-= Advanced reconstruction settings (for testing) =-=-\n", headerFont);
        gd.addCheckbox("Intensity weighting", prefs.get("intWeighting", true));
        gd.addCheckbox("Macro-pixel patterning correction", prefs.get("doMPmapCorrection", true));
        gd.addCheckbox("Show image stabilisation scatter plot", prefs.get("showImStabPlot", false));

        gd.addMessage("-=-= Advanced display settings =-=-\n", headerFont);
        gd.addCheckbox("Show gradients", prefs.get("showGradient", true));
        gd.addCheckbox("Show interpolated gradients", prefs.get("showIntGradient", true));

        gd.addHelp("https://www.youtube.com/watch?v=otCpCn0l4Wo"); // it's Hammer time

        MyDialogListenerAdvancedGUI dl = new MyDialogListenerAdvancedGUI(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            chosenDeviceName = prefs.get("chosenDeviceName", "Default device");
            maxMemoryGPU = prefs.get("maxMemoryGPU", 500);
            blockSize = (int) prefs.get("blockSize", 20000);
            intWeighting = prefs.get("intWeighting", true);
            doMPmapCorrection = prefs.get("doMPmapCorrection", true);
            showImStabPlot = prefs.get("showImStabPlot", false);
            showGradient = prefs.get("showGradient", true);
            showIntGradient = prefs.get("showIntGradient", true);

            // re-initialises to how it was before entering advanced GUI
            writeToDiskToUse = writeToDiskTemp;
            previousWriteToDisk = writeToDiskTemp;
        } else {
            writeToDiskToUse = writeToDiskTemp;
            prefs.set("chosenDeviceName", chosenDeviceName);
            prefs.set("maxMemoryGPU", maxMemoryGPU);
            prefs.set("blockSize", blockSize);
            prefs.set("intWeighting", intWeighting);
            prefs.set("doMPmapCorrection", doMPmapCorrection);
            prefs.set("showImStabPlot", showImStabPlot);
            prefs.set("showGradient", showGradient);
            prefs.set("showIntGradient", showIntGradient);
            prefs.save();
        }

    }

    // --- Advanced GUI Dialog listener ---
    class MyDialogListenerAdvancedGUI implements DialogListener {
        @Override
        public boolean dialogItemChanged(GenericDialog gd, AWTEvent awtEvent) {
            return grabSettingsAdvancedGUI(gd);
        }
    }


    // --- Grab settings from main GUI ---
    private boolean grabSettingsAdvancedGUI(GenericDialog gd) {

        chosenDeviceName = gd.getNextChoice();
        maxMemoryGPU = (float) gd.getNextNumber();
        blockSize = (int) gd.getNextNumber();

        boolean goodToGo = calculatenFrameOnGPU();

        writeToDiskTemp = gd.getNextBoolean();
        intWeighting = gd.getNextBoolean();
        doMPmapCorrection = gd.getNextBoolean();
        showImStabPlot = gd.getNextBoolean();
        showGradient = gd.getNextBoolean();
        showIntGradient = gd.getNextBoolean();

        if (writeToDiskTemp && !previousWriteToDisk) {
            pathToDisk = IJ.getDirectory("");
            if (pathToDisk == null) writeToDiskTemp = false;
        }

        previousWriteToDisk = writeToDiskTemp;
        return goodToGo;
    }


    // --- Predict memory usage in MB ---
    private double[] predictMemoryUsed(int nFrameOnGPU) { // TODO: adjust this for 3D case (more memory needed !)

        double[] memUsed = new double[2];

        // Memory on GPU ---- (estimated roughly, especially in teh case of 3D)
        memUsed[0] = 0;
        memUsed[0] += width * height * nFrameOnGPU; // clBufferPx
        memUsed[0] += nFrameOnGPU; // clBufferShiftX
        memUsed[0] += nFrameOnGPU; // clBufferShiftY
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGx
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGy
        if (do3DSRRF) memUsed[0] += width * height * nFrameOnGPU; // clBufferGz
        memUsed[0] += gradMag * gradMag * width * height * nFrameOnGPU; // clBufferGxInt
        memUsed[0] += gradMag * gradMag * width * height * nFrameOnGPU; // clBufferGyInt
        if (do3DSRRF) memUsed[0] += gradMag * gradMag * width * height * nFrameOnGPU; // clBufferGzInt

        memUsed[0] += nRecons * width * height * magnification * magnification; // clBufferRGC
        memUsed[0] += width * height * magnification * magnification; // clBufferInt

        // Memory on CPU ---- (largely underestimated)
        memUsed[1] = 0;
        memUsed[1] += width * height * nSlices; // RAM raw data
        if (do3DSRRF){
            memUsed[1] += nRecons * nSRRFframe * width * height * magnification * magnification * magnification; // clBufferSRRF
            memUsed[1] += nSRRFframe * width * height * magnification * magnification * magnification; // clBufferInt
        }
        else{
            memUsed[1] += nRecons * nSRRFframe * width * height * magnification * magnification; // clBufferSRRF
            memUsed[1] += nSRRFframe * width * height * magnification * magnification; // clBufferInt
        }

        memUsed[0] *= Float.SIZE / 8000000d; // convert to MB
        memUsed[1] *= Float.SIZE / 8000000d; // TODO: make a difference between Mac and PC?


        return memUsed;
    }

    // -- Save user's last parameter set --
    private void savePreferences() {

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("sensitivity", sensitivity);
        prefs.set("nFrameForSRRFfromUser", nFrameForSRRFfromUser);
        prefs.set("correctVibration", correctVibration);

        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateVAR", calculateVAR);
        prefs.set("calculateTAC2", calculateTAC2);
        prefs.set("getInterpolatedImage", getInterpolatedImage);

        prefs.set("doRollingAnalysis", doRollingAnalysis);
        prefs.set("frameGapFromUser", frameGapFromUser);

        prefs.set("chosenDeviceName", chosenDeviceName);
        prefs.set("maxMemoryGPU", maxMemoryGPU);
        prefs.set("blockSize", blockSize);

        prefs.set("showGradient", showGradient);
        prefs.set("showIntGradient", showIntGradient);

        prefs.set("do3DSRRF", do3DSRRF);

        prefs.save();

    }

    // --- Check the settings ---
    private boolean calculatenFrameOnGPU() {

        boolean goodToGo = false;

        int maxnFrameOnGPU = 0;
        double[] memUsed = new double[2];
        while (memUsed[0] < (double) maxMemoryGPU) {
            maxnFrameOnGPU = maxnFrameOnGPU + 1;
            memUsed = predictMemoryUsed(maxnFrameOnGPU);
        }

        maxnFrameOnGPU = maxnFrameOnGPU - 1;

        if (maxnFrameOnGPU > 0) {
            goodToGo = true;
            nFrameOnGPU = (int) math.ceil((float) nFrameForSRRFtoUse / (float) maxnFrameOnGPU);
            nFrameOnGPU = (int) math.ceil(((float) nFrameForSRRFtoUse / (float) nFrameOnGPU));
            nFrameOnGPU = min(nFrameForSRRFtoUse, nFrameOnGPU);
            memUsed = predictMemoryUsed(nFrameOnGPU);

            IJ.showStatus("eSRRF - Number of frames on device: " + (nFrameOnGPU) + " (" + (float) Math.ceil(memUsed[0]*100)/100 + " MB)");
        } else {
            memUsed = predictMemoryUsed(1);
            IJ.showStatus("eSRRF - Minimum device memory: " + (float) Math.ceil(memUsed[0]*100) / 100 + "MB");
        }

        // Check for RAM on computer
        if ((float) memUsed[1] > maxMemoryRAMij) {
            goodToGo = false;
            IJ.showStatus("eSRRF - Max RAM available exceeded: (" + ( (float) Math.ceil(memUsed[1]*100) / 100) + "MB vs. " + ((float) Math.round(maxMemoryRAMij*100) / 100) + "MB)");
            if (!writeSuggestOKeyed) IJ.showMessage("Results will likely exceed current RAM capacity. Consider increasing RAM for ImageJ or Write to disk (Advanced Settings) or else...");
            writeSuggestOKeyed = true;
        }

        return goodToGo;
    }


    // -- Convert time to string --
    private String timeToString(double time){

        String timeString;
        int _h = (int) (time / 3600);
        int _m = (int) (((time % 86400) % 3600) / 60);
        int _s = (int) (((time % 86400) % 3600) % 60);
        if (_h > 0) timeString = _h+"h "+_m+"m "+_s+"s";
        else {
            if (_m > 0) timeString = _m+"m "+_s+"s";
            else timeString = _s+"s";
        }

        return timeString;

    }

}
