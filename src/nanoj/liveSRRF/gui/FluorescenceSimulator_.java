package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.Binner;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core2.NanoJPrefs;

import java.util.LinkedHashMap;
import java.util.Map;

import static nanoj.core.java.array.ArrayMath.sum;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.core2.NanoJRandomNoise.gaussianValue;
import static nanoj.core2.NanoJRandomNoise.poissonValue;

import java.util.Random;



public class FluorescenceSimulator_ implements PlugIn {

    private ImagePlus imp;
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());


    public void run(String arg) {

        // Get raw data
        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        if ((imp.getBitDepth() == 24) || (imp.getBitDepth() == 32)){
            IJ.log("Please provide a 8-bit or 16-bit image as ground truth.");
            return;
        }

        Calibration cal = imp.getCalibration();
        if (cal.pixelHeight != cal.pixelWidth) return;
        float simPixelSize = (float) cal.pixelHeight;
        if ((cal.getUnit().equals("micron")) || (cal.getUnit().equals("um"))) simPixelSize *= 1000; // convert to nm // TODO: this doesnt work for some reason
        IJ.log("Ground truth pixel size: "+simPixelSize+" nm");

        ImageStack imsGT = imp.getStack();
        int[][] xy = getMoleculePositionsFromGroundTruth(imsGT);

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Fluorescence simulator");
        gd.addMessage("Microscope parameters");
        gd.addNumericField("Emission wavelength (nm):", prefs.get("lambda", 700), 2); // TODO: save prefs
        gd.addNumericField("NA:", prefs.get("NA", 1.15f), 2);
        gd.addNumericField("Number of frames: ", prefs.get("nFrames", 100), 0);
        gd.addNumericField("Exposure time (ms): ", prefs.get("exposure",10.0f), 1);
        gd.addNumericField("Pixel size (nm): ", prefs.get("pixelSize",100), 1);

        gd.addMessage("Camera parameters");
        gd.addNumericField("Gain (ADC/e-): ", prefs.get("cameraGain",2), 2);
        gd.addNumericField("Noise St.Dev. (e-): ", prefs.get("stdBgNoise", 2), 2);
        gd.addNumericField("Baseline level (ADC): ", prefs.get("cameraBaseline",100), 1);

        gd.addMessage("Fluorophore parameters");
        gd.addNumericField("kON (s^-1):", prefs.get("kON",10), 2);
        gd.addNumericField("kOFF (s^-1):", prefs.get("kOFF",10), 2);
        gd.addNumericField("Photon flux (s^-1):", prefs.get("photonFlux", 1000), 1);

        gd.addMessage("Simulation parameters");
        gd.addNumericField("Time division (max. 100): ", prefs.get("timeDivision",20), 0);
        gd.addCheckbox("Show background stack", prefs.get("showBgStack",false));
        gd.addCheckbox("Show simulated stack", prefs.get("showSimStack",false));
        gd.addCheckbox("Disable Poisson noise", prefs.get("disablePoissonNoise",false));
        gd.addCheckbox("Show individual traces (!)", prefs.get("showTraces", false));

        gd.showDialog();

        float lambda = (float) gd.getNextNumber();
        float NA = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        double exposure = gd.getNextNumber(); // exposure needs to be double as it's divided by 1000 and gets messed up as float
        float pixelSize = (float) gd.getNextNumber();

        float cameraGain = (float) gd.getNextNumber();
        float stdBgNoise = (float) gd.getNextNumber();
        float cameraBaseline = (float) gd.getNextNumber();

        float kON = (float) gd.getNextNumber();
        float kOFF = (float) gd.getNextNumber();
        float photonFlux = (float) gd.getNextNumber();

        int timeDivision = (int) gd.getNextNumber();
        boolean showBgStack = gd.getNextBoolean();
        boolean showSimStack = gd.getNextBoolean();
        boolean disablePoissonNoise = gd.getNextBoolean();
        boolean showTraces = gd.getNextBoolean();

        prefs.set("lambda", (float) lambda);
        prefs.set("NA", (float) NA);
        prefs.set("nFrames", (int) nFrames);
        prefs.set("exposure", (float) (exposure));
        prefs.set("pixelSize", (float) pixelSize);
        prefs.set("cameraGain", (float) cameraGain);
        prefs.set("stdBgNoise", (float) stdBgNoise);
        prefs.set("cameraBaseline", (float) cameraBaseline);
        prefs.set("kON", (float) kON);
        prefs.set("kOFF", (float) kOFF);
        prefs.set("photonFlux", (float) photonFlux);
        prefs.set("timeDivision", timeDivision);
        prefs.set("showBgStack", showBgStack);
        prefs.set("showSimStack", showSimStack);
        prefs.set("disablePoissonNoise", disablePoissonNoise);
        prefs.set("showTraces", showTraces);

        int wSim = imp.getWidth();
        int hSim = imp.getHeight();
        int nMolecules = xy.length;
        exposure /= 1000.0f; // convert to seconds

        double simDuration = (double) nFrames * exposure; // in seconds
        double sigmaPSF = 0.21*lambda/(NA*simPixelSize); // in simulated pixelsize
        double photonsPerTimeDivision = photonFlux*exposure/timeDivision;
        int shrinkFactor = (int) (pixelSize/simPixelSize);

        // Check if there's enough memory
        float maxMemoryRAMij = (float) IJ.maxMemory()/1e6f;  // in MB
        float neededRAM = (float) nMolecules*nFrames /1e6f; // in MB: intensityTraces
        neededRAM += (float) wSim*hSim*4 /1e6f; // in MB: Raw ground truth data
        neededRAM += (float) wSim*hSim*4 /1e6f; // in MB: Simulated single frame at ground truth scale
        neededRAM += (float) wSim*hSim*4 /1e6f; // in MB: Duplicated simulated single frame at ground truth scale (prior to binning) ?? not sure whether this is actually taking space
        neededRAM += (float) nFrames * (float) wSim/shrinkFactor * (float) hSim/shrinkFactor /1e6f; // in MB: xy coordinates from Ground truth
        neededRAM += (float) nMolecules*2 /1e6f; // in MB: xy coordinates from Ground truth

        IJ.log("-------------------");
        if (neededRAM < 1000) IJ.log("RAM required: "+neededRAM+" MB");
        else IJ.log("RAM required: "+(neededRAM/1000)+" GB");

        if (neededRAM > 0.9*maxMemoryRAMij) {
            IJ.log("ABORTED: not enough RAM to compute.");
            return;
        }
        IJ.log("-------------------");


        IJ.log("kON: "+kON+" s^-1 (Tau_OFF: "+1000*(1/kON)+"ms)");
        IJ.log("kOFF: "+kOFF+" s^-1 (Tau_ON: "+1000*(1/kOFF)+"ms)");
        IJ.log("Photon flux: "+photonFlux+" photons/s");
        IJ.log("Acquisition time: "+simDuration+"s ("+nFrames+" frames at "+exposure*1000+"ms exposure)");
        IJ.log("PSF sigma: "+(sigmaPSF*simPixelSize)+" nm");

        IJ.log("Creating variable for single molecule traces...");
        byte[] thisTrace;
        byte[][] intensityTraces = new byte[nMolecules][nFrames]; // minimising memory consumption here: this is the big boy
        IJ.log("RAM: " + IJ.freeMemory());

        IJ.log("Simulating time traces...");
        for (int i = 0; i < nMolecules; i++) {
            IJ.showProgress(i, nMolecules);
            thisTrace = getTimeTraceSingleMolecule(kON, kOFF, timeDivision, exposure, nFrames, showTraces);
            for (int j = 0; j < nFrames; j++) {
                intensityTraces[i][j] = thisTrace[j];
            }
        }

        ImageStack imsSim = new ImageStack();
        ImageStack imsFullSim = new ImageStack();
        ImageStack imsBg = new ImageStack();

        FloatProcessor fpSim;
        FloatProcessor fpSimBinned, fpBg;
        float[] pixels;
        Binner binner = new Binner();
        // Run garbage collector
        System.gc();

        IJ.log("Generating fluorescence images...");

        for (int f = 0; f < nFrames; f++) {

            // Run garbage collector
            System.gc();
            IJ.showProgress(f, nFrames);
            pixels = new float[wSim*hSim];

            for (int m = 0; m < nMolecules; m++) {
                pixels[xy[m][0]*wSim + xy[m][1]] += (float) intensityTraces[m][f] * photonsPerTimeDivision; // in photons
            }

            fpSim = new FloatProcessor(wSim, hSim, pixels);
            fpSim.blurGaussian(sigmaPSF);
            if (showSimStack) imsSim.addSlice(fpSim);

//            ipBinned = fpSim.duplicate(); // this also converts back to ImageProcessor, necessary for the Binner to work
//            ipBinned = binner.shrink(ipBinned, shrinkFactor, shrinkFactor, Binner.SUM);
//            fpSimBinned = ipBinned.duplicate().convertToFloatProcessor();

            fpSimBinned = binner.shrink(fpSim, shrinkFactor, shrinkFactor, Binner.SUM).convertToFloatProcessor();
//            fpSimBinned = fpSim.resize(w, h).convertToFloatProcessor();  // this doesn't do what we want. Doesn't add but decimates
            if (!disablePoissonNoise) fpSimBinned = addPoissonNoise(fpSimBinned.convertToFloatProcessor()).convertToFloatProcessor();
            fpSimBinned.multiply(cameraGain);

            fpBg = generateBackgroundGaussianNoise(fpSimBinned.getWidth(), fpSimBinned.getHeight(), cameraBaseline, stdBgNoise*cameraGain);
            if (showBgStack) imsBg.addSlice(fpBg);
            imsFullSim.addSlice(add(fpSimBinned.convertToFloatProcessor(), fpBg.convertToFloatProcessor()));
        }

        if (showSimStack) {
            ImagePlus impSim = new ImagePlus("Simulated stack", imsSim);
            Calibration calSim = impSim.getCalibration();
            calSim.setUnit("um");
            calSim.pixelHeight = simPixelSize/1000;
            calSim.pixelWidth  = simPixelSize/1000;
            impSim.show();
        }

        if (showBgStack) {
            ImagePlus impBg = new ImagePlus("Background stack", imsBg);
            Calibration calBg = impBg.getCalibration();
            calBg.setUnit("um");
            calBg.pixelHeight = pixelSize / 1000;
            calBg.pixelWidth = pixelSize / 1000;
            impBg.show();
        }

        ImagePlus impSimBinned = new ImagePlus(imp.getShortTitle()+" - Fluorescence stack", imsFullSim);
        Calibration calStack = impSimBinned.getCalibration();
        calStack.setUnit("um");
        calStack.pixelHeight = pixelSize/1000;
        calStack.pixelWidth  = pixelSize/1000;
        impSimBinned.show();

        // Run garbage collector
        System.gc();

        IJ.log("---------------------");
        IJ.log("All done.");

    }

    public final int[][] getMoleculePositionsFromGroundTruth(ImageStack imsGT){

        IJ.log("Getting molecule positions from the ground truth image...");
        float[] pixels = (float[]) imsGT.getProcessor(1).convertToFloatProcessor().getPixels();
        int nMolecules = 0;
        for (int i = 0; i < pixels.length; i++) {
            nMolecules += (int) pixels[i];
        }
        IJ.log("Number of molecules: "+nMolecules);
        int[][] xy = new int[nMolecules][2];
        double[] xArray = new double[nMolecules];
        double[] yArray = new double[nMolecules];

        int moleculeID = 0;
        for (int xi = 0; xi < imsGT.getWidth(); xi++) {
            for (int yi = 0; yi < imsGT.getHeight(); yi++) {
                int nMoleculeInPixel = (int) pixels[xi*imsGT.getWidth() + yi];
                for (int i = 0; i < nMoleculeInPixel; i++) {
                    xy[moleculeID][0] = xi;
                    xy[moleculeID][1] = yi;
                    xArray[moleculeID] = xi;
                    yArray[moleculeID] = yi;
                    moleculeID++;
                }
            }
        }

        // ---- Creating the NanoJ table and save it ----
        Map<String, double[]> data = new LinkedHashMap<>();
        data.put("X-position (pixels)", xArray);
        data.put("Y-position (pixels)", yArray);

        ResultsTable rt = dataMapToResultsTable(data);
        rt.show("Calibration-Table");

//        try {
//            saveNanoJTable(filePath+" - MFM_CalibrationTable.njt", rt);
//        } catch (IOException e) {
//            e.printStackTrace();
//        }

        return xy;
    }


    public static byte[] getTimeTraceSingleMolecule(float kON, float kOFF, int timeDivision, double exposure, int nFrames, boolean showTraces){

//        // ---- user set -----
//        boolean showTraces = true;
//        // -------------------

        double simTime = exposure/timeDivision;
        double duration = nFrames*exposure;

        Random rnd = new Random();
        boolean ON;
        double p = rnd.nextDouble();

        // Set initial probability of being ON at the start
        ON = (p < kON/(kON + kOFF));
        int traceLength = (int) (duration/simTime);
//        IJ.log("Tracelength: "+traceLength);

        byte[] switchingTrace = new byte[traceLength];
        double[] switchingTraceDouble = null;
        double[] intensityTraceDouble = null;
        if (showTraces) {
            switchingTraceDouble = new double[traceLength];
            intensityTraceDouble = new double[nFrames];
        }

        int nextTimeChange;
        int traceID = 0;
        int nSwitches = 0;

        // Generate the temporally fine trace of switching
        while ((traceID < traceLength)) {
//            IJ.log("Current trace ID: "+traceID);
//            IJ.log("Molecule is "+ON);
            p = rnd.nextDouble();
            if (ON) nextTimeChange = (int) (-Math.log(p) / (kOFF * simTime)); // in unit of simTime
            else nextTimeChange = (int) (-Math.log(p) / (kON * simTime)); // in unit of simTime
//            IJ.log("Next time change: "+(nextTimeChange*simTime)+" s");

            int i = 0;
            while ((i < nextTimeChange) && (traceID < traceLength)) {
                if (ON) {
                    switchingTrace[traceID] = 1;
                    if (showTraces) switchingTraceDouble[traceID] = 1;
                }
                i++;
                traceID++;
            }
            ON = !ON;
            nSwitches ++;
        }

        // Bin the time trace to the frame size already (huge saving in memory, a factor of timeDivision)
        byte[] intensityTrace = new byte[nFrames];
        for (int f = 0; f < nFrames; f++) {
            for (int i = 0; i < timeDivision; i++) {
                intensityTrace[f] += switchingTrace[f*timeDivision + i]; // byte should not exceed 127!!! TODO: adjust to use full 8-bits by starting at -128?
                if (showTraces) intensityTraceDouble[f] += (double) switchingTrace[f*timeDivision + i];
            }
        }

        // Plot the switching trace if necessary
        if (showTraces) {
            Plot blinkingTracePlot = new Plot("Blink trace plot", "Time (time division)", "Switching state");
            blinkingTracePlot.add("line", switchingTraceDouble);
            blinkingTracePlot.show();

            Plot IntensityTracePlot = new Plot("Intensity trace plot", "Time (frames)", "Intensity state");
            IntensityTracePlot.add("line", intensityTraceDouble);
            IntensityTracePlot.show();

            // Calculating the ON and OFF times from the trace and print it
            double avONtime = (float) sum(switchingTraceDouble)*simTime*1000/(nSwitches/2); // in ms
            double avOFFtime = 1000*(duration - (float) sum(switchingTraceDouble)*simTime)/(nSwitches/2);  // in ms
            IJ.log("Number of blinks: "+(nSwitches/2));
            IJ.log("Average ON time: "+avONtime+" ms");
            IJ.log("Average OFF time: "+avOFFtime+" ms");
        }

        return intensityTrace;
    }


    public static FloatProcessor addPoissonNoise(FloatProcessor fp){

        float[] pixels = (float[]) fp.getPixels();
        float[] pixelsPoisson = new float[pixels.length];
        for (int i = 0; i < pixels.length; i++) pixelsPoisson[i] = (float) poissonValue(pixels[i]);
        return new FloatProcessor(fp.getWidth(), fp.getHeight(), pixelsPoisson);
    }

    public static FloatProcessor generateBackgroundGaussianNoise(int width, int height, float mean, float std){

        float[] pixelsGaussian = new float[width*height];
        for (int i = 0; i < width*height; i++) pixelsGaussian[i] = (float) gaussianValue(mean, std);
        return new FloatProcessor(width, height, pixelsGaussian);
    }

}
