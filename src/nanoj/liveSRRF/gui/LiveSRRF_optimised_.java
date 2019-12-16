package nanoj.liveSRRF.gui;

import com.jogamp.opencl.CLDevice;
import ij.*;
import ij.gui.DialogListener;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import nanoj.core.java.io.zip.SaveFileInZip;
import nanoj.core.java.io.zip.virtualStacks.FullFramesVirtualStack;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.core2.NanoJUsageTracker;
import nanoj.liveSRRF.XYShiftCalculator;
import org.python.modules.math;
import nanoj.liveSRRF.liveSRRF_CL;

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
            width,
            height,
            blockSize,
            nGPUloadPerSRRFframe;

    private float fwhm, maxMemoryGPU,
            maxMemoryRAMij = (float) IJ.maxMemory()/1e6f; // maximum RAM set for Fiji in MB

    private boolean correctVibration,
            calculateAVG,
            calculateSTD,
            getInterpolatedImage,
            writeToDiskToUse,
            writeToDiskTemp,
            doRollingAnalysis,
            previousWriteToDisk = false,
            previousAdvSettings = false,
            writeSuggestOKeyed = false,
            doMPmapCorrection,
            showImStabPlot,
            intWeighting;

    private final String LiveSRRFVersion = "v1.10";
    private String pathToDisk = "",
            fileName,
            chosenDeviceName,
            chosenTemporalAnalysis;

    private String[] deviceNames,
            temporalAnalysis = {"AVG","STD","Both AVG and STD"};

    private float[] shiftX, shiftY;

    // Image formats
    private ImagePlus imp;

    private ImagePlus impSRRFavg,
            impSRRFstd,
            impRawInterpolated;


    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();
    private liveSRRF_CL liveSRRF;
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
        IJ.log("liveSRRF " + LiveSRRFVersion); //TODO: what's taking so long on Wolverine?
        IJ.log("Max RAM available: "+ (float) Math.round(maxMemoryRAMij*100)/100 + " MB");

        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));

        // Initialize the liveSRRF class and check the devices
        liveSRRF = new liveSRRF_CL();
        liveSRRF.checkDevices();
        CLDevice[] allDevices = liveSRRF_CL.allCLdevices;

        // Initializing string for device choice
        deviceNames = new String[allDevices.length + 1];
        deviceNames[0] = "Default device";

        for (int i = 1; i <= allDevices.length; i++) {
            deviceNames[i] = allDevices[i - 1].getName();
        }

        // Main GUI ----
        boolean cancelled = mainGUI();
        if (cancelled) {
            liveSRRF.release();
            IJ.log("Cancelled by user.");
            return;
        }

        // Check if something has gone wrong with the momory
        if (!calculatenFrameOnGPU()){
            liveSRRF.release();
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

        liveSRRF.initialise(width, height, magnification, fwhm, sensitivity, nFrameOnGPU, nFrameForSRRFtoUse, blockSize, chosenDevice, intWeighting, doMPmapCorrection);

        shiftX = new float[nFrameForSRRFtoUse];
        shiftY = new float[nFrameForSRRFtoUse];

        XYShiftCalculator shiftCalculator = new XYShiftCalculator(imp, prof);
        ImageStack imsBuffer;

        ImageStack imsRawData;
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);
        ImageStack imsRawInterpolated = new ImageStack(width * magnification, height * magnification);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

        IJ.log("Distance units: "+cal.getUnit());

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

            liveSRRF.loadShiftXYGPUbuffer(shiftX, shiftY);

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

                userPressedEscape = liveSRRF.calculateSRRF(imsRawData); // resets the local GPU load frame counter

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
                        impTemp = new ImagePlus("\"liveSRRF (AVG) - frame=" + r + ".tif", imsBuffer.getProcessor(1));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("liveSRRF (AVG) - frame=" + r, impTemp);
                    }
                    if (calculateSTD) {
                        impTemp = new ImagePlus("\"liveSRRF (STD) - frame=" + r + ".tif", imsBuffer.getProcessor(2));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("liveSRRF (STD) - frame=" + r, impTemp);
                    }
                    if (getInterpolatedImage) {
                        impTemp = new ImagePlus("\"liveSRRF (INT) - frame=" + r + ".tif", imsBuffer.getProcessor(3));
                        impTemp.setCalibration(cal);
                        saveFileInZip.addTiffImage("liveSRRF (INT) - frame=" + r, impTemp);
                    }
                } catch (IOException e) {
                    IJ.error("Whoops, it seems that there was a problem with saving to disk... (insert sad face here).");
                    e.printStackTrace();
                }


            } else {
                if (calculateAVG) imsSRRFavg.addSlice(imsBuffer.getProcessor(1));
                if (calculateSTD) imsSRRFstd.addSlice(imsBuffer.getProcessor(2));
                if (getInterpolatedImage) imsRawInterpolated.addSlice(imsBuffer.getProcessor(3));
            }
            IJ.log("RAM used: " + IJ.freeMemory());
        }
        // End looping trough SRRF frames --------------------------------------------


        // Release the GPU
        liveSRRF.release();

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

                    impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG)", vsimsSRRFavg);
                    impSRRFavg.setCalibration(cal);
                    IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                    impSRRFavg.show();
                }

                if (calculateSTD) {
                    FullFramesVirtualStack vsimsSRRFstd = new FullFramesVirtualStack(fileName, true);
                    for (int f = vsimsSRRFstd.getSize(); f > 0; f--) {
                        if (!vsimsSRRFstd.getSliceLabel(f).contains("STD")) vsimsSRRFstd.deleteSlice(f);
                    }

                    impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD)", vsimsSRRFstd);
                    impSRRFstd.setCalibration(cal);
                    IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
                    impSRRFstd.show();
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
                impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG)", imsSRRFavg);
                impSRRFavg.setCalibration(cal);
                IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
                impSRRFavg.show();
            }

            if (calculateSTD) {
                impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD)", imsSRRFstd);
                impSRRFstd.setCalibration(cal);
                IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
                impSRRFstd.show();
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


//        IJ.log("-------------------------------------");
//        IJ.log(prof.report());
    }


    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    // -- Main GUI --
    private boolean mainGUI() {
        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("liveSRRF " + LiveSRRFVersion);
        gd.addMessage("-=-= SRRF parameters =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);
        gd.addNumericField("Radius (pixels, default: 1.5)", prefs.get("fwhm", (float) 1.5), 2);
        gd.addNumericField("Sensitivity (default: 1)", prefs.get("sensitivity", 1), 0);
        gd.addNumericField("# frames for SRRF (0 = auto)", prefs.get("nFrameForSRRFfromUser", 0), 0);
        gd.addCheckbox("Vibration correction", prefs.get("correctVibration", false));

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addChoice("Temporal analysis", temporalAnalysis, prefs.get("chosenTemporalAnalysis", temporalAnalysis[2]));
//        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
//        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("calculateSTD", false));
        gd.addCheckbox("Wide-field interpolation (default: off)", prefs.get("getInterpolatedImage", false));

        gd.addMessage("-=-= Rolling analysis =-=-\n", headerFont);
        gd.addCheckbox("Perform rolling analysis (default: off)", prefs.get("doRollingAnalysis", false));
        gd.addNumericField("# frame gap between SR frame (0 = auto)", prefs.get("frameGapFromUser", 0), 0);

        gd.addMessage("-=-= Advanced settings =-=-\n", headerFont);
        gd.addCheckbox("Show advanced settings", false);

        gd.addHelp("https://www.youtube.com/watch?v=PJQVlVHsFF8"); // If you're hooked on a feeling

        MyDialogListenerMainGUI dl = new MyDialogListenerMainGUI(); // this serves to estimate a few indicators such as RAM usage
        gd.addDialogListener(dl);

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

        chosenTemporalAnalysis = gd.getNextChoice();
        if (chosenTemporalAnalysis.equals(temporalAnalysis[0])){
            calculateAVG = true;
            calculateSTD = false;
        }
        else if(chosenTemporalAnalysis.equals(temporalAnalysis[1])){
            calculateAVG = false;
            calculateSTD = true;
        }
        else{
            calculateAVG = true;
            calculateSTD = true;
        }

//        calculateAVG = gd.getNextBoolean();
//        calculateSTD = gd.getNextBoolean();
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

        nGPUloadPerSRRFframe = (int) math.ceil((float) nFrameForSRRFtoUse / (float) nFrameOnGPU);

        return goodToGo;
    }


    // -- Advanced settings GUI --
    private void advancedSettingsGUI() {

        double[] memUsed = predictMemoryUsed(1);

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        GenericDialog gd = new GenericDialog("liveSRRF - Advanced settings");

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

        // Memory on GPU ----
        memUsed[0] = 0;
        memUsed[0] += width * height * nFrameOnGPU; // clBufferPx
        memUsed[0] += nFrameOnGPU; // clBufferShiftX
        memUsed[0] += nFrameOnGPU; // clBufferShiftY
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGx
        memUsed[0] += width * height * nFrameOnGPU; // clBufferGy
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGxInt
        memUsed[0] += 4 * width * height * nFrameOnGPU; // clBufferGyInt
        memUsed[0] += 2 * width * height * magnification * magnification; // clBufferRGC
        memUsed[0] += width * height * magnification * magnification; // clBufferInt

        // Memory on CPU ---- (largely underestimated)
        memUsed[1] = 0;
        memUsed[1] += width * height * nSlices; // clBufferPx
        memUsed[1] += 2 * nSRRFframe * width * height * magnification * magnification; // clBufferSRRF
        memUsed[1] += nSRRFframe * width * height * magnification * magnification; // clBufferInt

        memUsed[0] *= Float.SIZE / 8000000d;
        memUsed[1] *= Float.SIZE / 8000000d;


        return memUsed;
    }

    // -- Save user's last parameter set --
    private void savePreferences() {

        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm", fwhm);
        prefs.set("nFrameForSRRFfromUser", nFrameForSRRFfromUser);
        prefs.set("sensitivity", sensitivity);
        prefs.set("correctVibration", correctVibration);

        prefs.set("chosenTemporalAnalysis", chosenTemporalAnalysis);
//        prefs.set("calculateAVG", calculateAVG);
//        prefs.set("calculateSTD", calculateSTD);
        prefs.set("getInterpolatedImage", getInterpolatedImage);

        prefs.set("doRollingAnalysis", doRollingAnalysis);
        prefs.set("frameGapFromUser", frameGapFromUser);

        prefs.set("chosenDeviceName", chosenDeviceName);
        prefs.set("maxMemoryGPU", maxMemoryGPU);
        prefs.set("blockSize", blockSize);

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

            IJ.showStatus("liveSRRF - Number of frames on device: " + (nFrameOnGPU) + " (" + (float) Math.ceil(memUsed[0]*100)/100 + " MB)");
        } else {
            memUsed = predictMemoryUsed(1);
            IJ.showStatus("liveSRRF - Minimum device memory: " + (float) Math.ceil(memUsed[0]*100) / 100 + "MB");
        }

        // Check for RAM on computer
        if ((float) memUsed[1] > maxMemoryRAMij) {
            goodToGo = false;
            IJ.showStatus("liveSRRF - Max RAM available exceeded: (" + ( (float) Math.ceil(memUsed[1]*100) / 100) + "MB vs. " + ((float) Math.round(maxMemoryRAMij*100) / 100) + "MB)");
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
