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
            nFramesRawData,
            width,
            height,
            nRecons,
            blockSize,
            nGPUloadPerSRRFframe;

    private float fwhm, maxMemoryGPU,
            maxMemoryRAMij = (float) IJ.maxMemory()/1e6f; // maximum RAM set for Fiji in MB

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
            showGradients,
            intWeighting;

    private final String LiveSRRFVersion = "v1.2d-fhi.14";
    private String pathToDisk = "",
            fileName,
            chosenDeviceName,
//            chosenTemporalAnalysis,
            thisGradientChoice;

    private String[] deviceNames,
//            temporalAnalysis = {"AVG","STD","Both AVG and STD"},
            gradientChoice = {"RobX", "3pPlus", "3pX"};

    private float[][] driftXY;

    // Image formats
    private ImagePlus imp,
            impSRRFavg,
            impSRRFvar,
            impSRRFtac2,
            impRawInterpolated;


    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private LiveSRRF_CL liveSRRF;
    private SaveFileInZip saveFileInZip;
    private CLDevice chosenDevice = null;

    // Tracker
    private String user = "FijiUpdater";
    private String version = "20180809-" + user;
    private NanoJUsageTracker tracker = new NanoJUsageTracker("NanoJ-LiveSRRF", version, "UA-61590656-4");


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

        nFramesRawData = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();

        // initialization of advanced settings (in case advanced settings are not selected)
        chosenDeviceName = prefs.get("chosenDeviceName", "Default device");
        maxMemoryGPU = prefs.get("maxMemoryGPU", 500);
        blockSize = (int) prefs.get("blockSize", 20000);
        writeToDiskToUse = false;
        intWeighting = prefs.get("intWeighting", true);
        doMPmapCorrection = prefs.get("doMPmapCorrection", true);
        showGradients = prefs.get("showGradients", false);


//        // Close the log: //TODO: this seems to sometimes make a different whether OpenCL runs or not (observation from Nvidia 1050, not repeatable)
        // TODO: but this causes the GUI to pop up when using liveSRRF in a macro
//        if (IJ.getLog() != null) {
//            selectWindow("Log");
//            IJ.run("Close");
//        }

        IJ.log("\\Clear");  // Clear the log window

        // Check saved preferences
//        IJ.log("---- Check preferences saved ----");
//        IJ.log("Magnification: "+prefs.get("magnification", 100));
//        IJ.log("FWHM: "+prefs.get("fwhm", 100));
//        IJ.log("Sensitivity: "+prefs.get("sensitivity", 100));
//        IJ.log("nFrameForSRRFfromUser: "+prefs.get("nFrameForSRRFfromUser", 100));


        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("LiveSRRF " + LiveSRRFVersion);
        IJ.log("Max RAM available: "+ (float) Math.round(maxMemoryRAMij*100)/100 + " MB");

        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));
        IJ.log("Pixel size: "+pixelSize+" "+pixelSizeUnit);

        // Initialize the liveSRRF class and check the devices
        liveSRRF = new LiveSRRF_CL();
        liveSRRF.checkDevices();
        CLDevice[] allDevices = LiveSRRF_CL.allCLdevices;
        nRecons = liveSRRF.nReconstructions;

        // Initializing string for device choice
        deviceNames = new String[allDevices.length + 1];
        deviceNames[0] = "Default device";

        for (int i = 1; i <= allDevices.length; i++) {
            deviceNames[i] = allDevices[i - 1].getName();
        }

        // Main GUI ----
        boolean cancelled = mainGUI();
        if (cancelled) {
//            liveSRRF.release(); // no need te release contex prior to initilisation of liveSRRF
            IJ.log("Cancelled by user.");
            return;
        }

        // Check if something has gone wrong with the memory
        if (!calculatenFrameOnGPU()){
//            liveSRRF.release(); // no need te release contex prior to initilisation of liveSRRF
            IJ.log("Problems with memory. Check advanced settings or else...");
            return;
        }

        // Get chosen device
        if (chosenDeviceName == null) chosenDeviceName = "Default device";
        for (CLDevice thisDevice : allDevices) {
            if (chosenDeviceName.equals(thisDevice.getName())) chosenDevice = thisDevice;
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
            fileName = pathToDisk + imp.getTitle() + " - liveSRRF.zip";
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

        liveSRRF.initialise(width, height, magnification, fwhm, sensitivity, nFrameOnGPU, nFrameForSRRFtoUse, blockSize, chosenDevice, intWeighting, doMPmapCorrection, thisGradientChoice);

        driftXY = new float[nFrameForSRRFtoUse][2];
//        driftY = new float[nFrameForSRRFtoUse];

        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp);
        ImageStack imsBuffer;

        ImageStack imsRawData;
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFvar = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFtac2 = new ImageStack(width * magnification, height * magnification);
        ImageStack imsRawInterpolated = new ImageStack(width * magnification, height * magnification);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

//        IJ.log("Distance units: "+cal.getUnit());

        boolean userPressedEscape;
        long loopStart = System.nanoTime();

        // For measuring the remaining time
        double singleLoadTime;
        double remainingTime;
        int currentNLoad;

        // Start looping trough SRRF frames --------------------------------------------
        for (int r = 1; r <= nSRRFframe; r++) {
            liveSRRF.resetFramePosition();  // resets only the global SRRF counter

            IJ.log("--------");
            IJ.log("SRRF frame: " + r + "/" + nSRRFframe);
            indexStartSRRFframe = (r - 1) * frameGapToUse + 1;
            IJ.log("Stack index start: " + indexStartSRRFframe);

            if (correctVibration) {
                shiftCalculator.calculateShiftArray(indexStartSRRFframe, nFrameForSRRFtoUse);
                for (int i = 0; i < nFrameForSRRFtoUse; i++) {
                    driftXY[i][0] = shiftCalculator.shiftX[i];
                    driftXY[i][1] = shiftCalculator.shiftY[i];
                }


                if (showImStabPlot) {
                    // Scatter plot
                    Plot scatterPlot = new Plot("x/y scatter plot (spectrum LUT-coded) #" + r, "x (pixels)", "y (pixels)");
                    for (int i = 0; i < nFrameForSRRFtoUse; i++) {
                        float[] x_temp = new float[1];
                        float[] y_temp = new float[1];
                        x_temp[0] = driftXY[i][0];
                        y_temp[0] = driftXY[i][1];
                        scatterPlot.setColor(Color.getHSBColor(i / (float) nFramesRawData, 1f, 1f)); // this corresponds to the spectrum LUT
                        scatterPlot.addPoints(x_temp, y_temp, Plot.CROSS);
                    }
                    scatterPlot.show();
                    scatterPlot.setLimitsToFit(true);
                }

            }

            liveSRRF.loadDriftXYGPUbuffer(driftXY);

            IJ.showProgress(r - 1, nSRRFframe);

            for (int l = 0; l < nGPUloadPerSRRFframe; l++) {

//                IJ.log("----------------------------");
                nFrameToLoad = min(nFrameOnGPU, nFrameForSRRFtoUse - nFrameOnGPU * l);
//                IJ.log("Number of frames to load: " + nFrameToLoad);
//                IJ.log("Index start: " + (indexStartSRRFframe + nFrameOnGPU * l));

                imsRawData = new ImageStack(width, height);
                for (int f = 0; f < nFrameToLoad; f++) {
                    imsRawData.addSlice(imsAllRawData.getProcessor(indexStartSRRFframe + nFrameOnGPU * l + f));
                }

                liveSRRF.prepareDataSRRF(imsRawData);
                userPressedEscape = liveSRRF.calculateSRRF(); // resets the local GPU load frame counter

                // Check if user is cancelling calculation
                if (userPressedEscape) {
                    liveSRRF.release();
                    IJ.log("-------------------------------------");
                    IJ.log("Reconstruction aborted by user.");
                    return;
                }

                // Estimating the remaining time
                currentNLoad = nGPUloadPerSRRFframe*(r-1) + l+1;
                singleLoadTime = ((System.nanoTime() - loopStart) / (float) currentNLoad) / 1e9;
                remainingTime = singleLoadTime * (nSRRFframe * nGPUloadPerSRRFframe - currentNLoad);
                IJ.showStatus("LiveSRRF - Remaining time: " + timeToString(remainingTime));
            }

            liveSRRF.readSRRFbuffer();
            imsBuffer = liveSRRF.imsSRRF;

//            impMPmap = new ImagePlus("MP map",liveSRRF.getMPmap());
//            impMPmap.setCalibration(cal);
//            IJ.run(impMPmap, "Enhance Contrast", "saturated=0.5");
//            impMPmap.show();

            if (writeToDiskToUse) {
                try {
                    if (calculateAVG) {
                        impTemp = new ImagePlus("\"LiveSRRF (AVG) - frame=" + r + ".tif", imsBuffer.getProcessor(1));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("LiveSRRF (AVG) - frame=" + r, impTemp);
                    }
                    if (calculateVAR) {
                        impTemp = new ImagePlus("\"LiveSRRF (VAR) - frame=" + r + ".tif", imsBuffer.getProcessor(2));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("LiveSRRF (VAR) - frame=" + r, impTemp);
                    }
                    if (calculateTAC2) {
                        impTemp = new ImagePlus("\"LiveSRRF (TAC2) - frame=" + r + ".tif", imsBuffer.getProcessor(3));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("LiveSRRF (TAC2) - frame=" + r, impTemp);
                    }
                    if (getInterpolatedImage) {
                        impTemp = new ImagePlus("\"LiveSRRF (INT) - frame=" + r + ".tif", imsBuffer.getProcessor(4));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("LiveSRRF (INT) - frame=" + r, impTemp);
                    }
                } catch (IOException e) {
                    IJ.error("Whoops, it seems that there was a problem with saving to disk... (insert sad face here).");
                    e.printStackTrace();
                }


            } else {
                if (calculateAVG) imsSRRFavg.addSlice(imsBuffer.getProcessor(1));
                if (calculateVAR) imsSRRFvar.addSlice(imsBuffer.getProcessor(2));
                if (calculateTAC2) imsSRRFtac2.addSlice(imsBuffer.getProcessor(3));
                if (getInterpolatedImage) imsRawInterpolated.addSlice(imsBuffer.getProcessor(4));
            }
            IJ.log("RAM used: " + IJ.freeMemory());
        }
        // End looping trough SRRF frames --------------------------------------------

        if (showGradients) { // TODO: small bug leading to not remembering choices of showGradient
            // Read off the gradient (single frame, which one needs to be checked)
            IJ.log("Reading off and displaying gradients...");
            System.out.println("Reading off and displaying gradients...");
            ImageStack imsGradient = liveSRRF.readGradientBuffers();
            ImagePlus impGradients = new ImagePlus("Gradients", imsGradient);
            impGradients.show();
        }

        // Release the GPU
        liveSRRF.release();
        liveSRRF.checkProfiler();

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

                    impSRRFavg = new ImagePlus(imp.getTitle() + " - LiveSRRF (AVG)", vsimsSRRFavg);
                    impSRRFavg.setCalibration(cal);
                    IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                    impSRRFavg.show();
                }

                if (calculateVAR) {
                    FullFramesVirtualStack vsimsSRRFvar = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFvar.getSize(); f > 0; f--) {
                        if (!vsimsSRRFvar.getSliceLabel(f).contains("VAR")) vsimsSRRFvar.deleteSlice(f);
                    }

                    impSRRFvar = new ImagePlus(imp.getTitle() + " - LiveSRRF (VAR)", vsimsSRRFvar);
                    impSRRFvar.setCalibration(cal);
                    IJ.run(impSRRFvar, "Enhance Contrast", "saturated=0.5");
                    impSRRFvar.show();
                }

                if (calculateTAC2) {
                    FullFramesVirtualStack vsimsSRRFtac2 = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFtac2.getSize(); f > 0; f--) {
                        if (!vsimsSRRFtac2.getSliceLabel(f).contains("TAC2")) vsimsSRRFtac2.deleteSlice(f);
                    }

                    impSRRFtac2 = new ImagePlus(imp.getTitle() + " - LiveSRRF (TAC2)", vsimsSRRFtac2);
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
                impSRRFavg = new ImagePlus(imp.getTitle() + " - LiveSRRF (AVG)", imsSRRFavg);
                impSRRFavg.setCalibration(cal);
                IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                impSRRFavg.show();
            }

            if (calculateVAR) {
                impSRRFvar = new ImagePlus(imp.getTitle() + " - LiveSRRF (VAR)", imsSRRFvar);
                impSRRFvar.setCalibration(cal);
                IJ.run(impSRRFvar, "Enhance Contrast", "saturated=0.5");
                impSRRFvar.show();
            }

            if (calculateTAC2) {
                impSRRFtac2 = new ImagePlus(imp.getTitle() + " - LiveSRRF (TAC2)", imsSRRFtac2);
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

        // Check saved preferences
//        IJ.log("---- Check preferences saved ----");
//        IJ.log("Magnification: "+prefs.get("magnification", 100));
//        IJ.log("FWHM: "+prefs.get("fwhm", 100));
//        IJ.log("Sensitivity: "+prefs.get("sensitivity", 100));
//        IJ.log("nFrameForSRRFfromUser: "+prefs.get("nFrameForSRRFfromUser", 100));



    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    // -- Main GUI --
    private boolean mainGUI() {
        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("LiveSRRF " + LiveSRRFVersion);
        gd.addMessage("-=-= SRRF parameters =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("Radius (pixels, default: 1.5)", prefs.get("fwhm", (float) 1.5), 2);
        gd.addNumericField("Sensitivity (default: 1)", prefs.get("sensitivity", 1), 0);
        gd.addNumericField("# frames for SRRF (0 = auto)", prefs.get("nFrameForSRRFfromUser", 0), 0);
        gd.addCheckbox("Vibration correction", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
//        gd.addChoice("Temporal analysis", temporalAnalysis, prefs.get("chosenTemporalAnalysis", temporalAnalysis[2]));
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("VAR reconstruction (default: off)", prefs.get("calculateVAR", false));
        gd.addCheckbox("TAC2 reconstruction (default: off)", prefs.get("calculateTAC2", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("getInterpolatedImage", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addCheckbox("Perform rolling analysis (default: off)", prefs.get("doRollingAnalysis", false));
        gd.addNumericField("# frame gap between SR frame (0 = auto)", prefs.get("frameGapFromUser", 0), 0);

        gd.addMessage("-=-= Advanced settings =-=-\n", headerFont);
        gd.addCheckbox("Show advanced settings", false);
        gd.addChoice("Gradient type", gradientChoice, gradientChoice[0]);
//        gd.addChoice("Gradient mag.", new int[] {1,2}, 1);

        gd.addHelp("https://www.youtube.com/watch?v=PJQVlVHsFF8"); // If you're hooked on a feeling

        MyDialogListenerMainGUI dl = new MyDialogListenerMainGUI(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);
        grabSettingsMainGUI(gd); // grab and assign the default values to the global variables (in case the dialog listener is not activated, pressed OK directly)

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

//        chosenTemporalAnalysis = gd.getNextChoice();
//        if (chosenTemporalAnalysis.equals(temporalAnalysis[0])){
//            calculateAVG = true;
//            calculateVAR = false;
//        }
//        else if(chosenTemporalAnalysis.equals(temporalAnalysis[1])){
//            calculateAVG = false;
//            calculateVAR = true;
//        }
//        else{
//            calculateAVG = true;
//            calculateVAR = true;
//        }

        calculateAVG = gd.getNextBoolean();
        calculateVAR = gd.getNextBoolean();
        calculateTAC2 = gd.getNextBoolean();
        getInterpolatedImage = gd.getNextBoolean();

        doRollingAnalysis = gd.getNextBoolean();
        frameGapFromUser = (int) gd.getNextNumber();

        // Calculate the parameters based on user inputs
        if (nFrameForSRRFfromUser == 0) nFrameForSRRFtoUse = nFramesRawData;
        else nFrameForSRRFtoUse = min(nFramesRawData, nFrameForSRRFfromUser);

        if (frameGapFromUser == 0 || !doRollingAnalysis) frameGapToUse = nFrameForSRRFtoUse;
        else frameGapToUse = frameGapFromUser;

        nSRRFframe = (int) ((float) (nFramesRawData - nFrameForSRRFtoUse) / (float) frameGapToUse) + 1;

        boolean goodToGo = calculatenFrameOnGPU();

        boolean showAdvancedSettings = gd.getNextBoolean();
        if (showAdvancedSettings && !previousAdvSettings) advancedSettingsGUI();
        previousAdvSettings = showAdvancedSettings;

        thisGradientChoice = gd.getNextChoice();

        nGPUloadPerSRRFframe = (int) math.ceil((float) nFrameForSRRFtoUse / (float) nFrameOnGPU);

        return goodToGo;
    }


    // -- Advanced settings GUI --
    private void advancedSettingsGUI() {

        double[] memUsed = predictMemoryUsed(1);

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        GenericDialog gd = new GenericDialog("LiveSRRF - Advanced settings");

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
        gd.addCheckbox("Show gradients", prefs.get("showGradients", false));


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
            showGradients = prefs.get("showGradients", false);

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
            prefs.set("showGradients", showGradients);
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
        showGradients = gd.getNextBoolean();

        if (writeToDiskTemp && !previousWriteToDisk) {
            pathToDisk = IJ.getDirectory("");
            if (pathToDisk == null) writeToDiskTemp = false;
        }

        previousWriteToDisk = writeToDiskTemp;
        return goodToGo;
    }


    // --- Predict memory usage in MB ---
    private double[] predictMemoryUsed(int nFrameOnGPU) {

        double[] memUsed = new double[2];

        // Memory on GPU ---- TODO: update with the FHT interpolation memory requirements
        memUsed[0] = 0;
        memUsed[0] += width*magnification * height*magnification * nFrameOnGPU; // clBufferPx
        memUsed[0] += 2*nFrameOnGPU; // clBufferDriftXY
        memUsed[0] += width*magnification * height*magnification * nFrameOnGPU; // clBufferGx
        memUsed[0] += width*magnification * height*magnification * nFrameOnGPU; // clBufferGy

        memUsed[0] += width * height * magnification * magnification; // clBufferPreviousFrame
        memUsed[0] += (nRecons+1) * width * height * magnification * magnification; // clBufferOut

        // Memory on CPU ---- (largely underestimated)
        memUsed[1] = 0;
        memUsed[1] += width * height * nFramesRawData; // raw data
        memUsed[1] += (nRecons+1) * nSRRFframe * width * height * magnification * magnification; // Results

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

//        prefs.set("chosenTemporalAnalysis", chosenTemporalAnalysis);
        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateVAR", calculateVAR);
        prefs.set("calculateTAC2", calculateTAC2);
        prefs.set("getInterpolatedImage", getInterpolatedImage);

        prefs.set("doRollingAnalysis", doRollingAnalysis);
        prefs.set("frameGapFromUser", frameGapFromUser);

        prefs.set("chosenDeviceName", chosenDeviceName);
        prefs.set("maxMemoryGPU", maxMemoryGPU);
        prefs.set("blockSize", blockSize);
        prefs.set("showGradient", showGradients);

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

            IJ.showStatus("LiveSRRF - Number of frames on device: " + (nFrameOnGPU) + " (" + (float) Math.ceil(memUsed[0]*100)/100 + " MB)");
        } else {
            memUsed = predictMemoryUsed(1);
            IJ.showStatus("LiveSRRF - Minimum device memory: " + (float) Math.ceil(memUsed[0]*100) / 100 + "MB");
        }

        // Check for RAM on computer
        if ((float) memUsed[1] > maxMemoryRAMij) {
            goodToGo = false;
            IJ.showStatus("LiveSRRF - Max RAM available exceeded: (" + ( (float) Math.ceil(memUsed[1]*100) / 100) + "MB vs. " + ((float) Math.round(maxMemoryRAMij*100) / 100) + "MB)");
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
