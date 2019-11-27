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
import ij.process.ImageProcessor;
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
        gd.addNumericField("Exposure time (ms): ", prefs.get("exposure",10), 1);
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
        gd.addNumericField("Time division: ", prefs.get("timeDivision",20), 0);
        gd.addCheckbox("Show background stack", prefs.get("showBgStack",false));
        gd.addCheckbox("Show simulated stack", prefs.get("showSimStack",false));
        gd.addCheckbox("Disable Poisson noise", prefs.get("disablePoissonNoise",false));


        gd.showDialog();

        float lambda = (float) gd.getNextNumber();
        float NA = (float) gd.getNextNumber();
        int nFrames = (int) gd.getNextNumber();
        float exposure = ((float) gd.getNextNumber())/1000; // convert to seconds
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

        prefs.set("lambda", (float) lambda);
        prefs.set("NA", (float) NA);
        prefs.set("nFrames", nFrames);
        prefs.set("exposure", (float) (1000*exposure));
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


        int nMolecules = xy.length;
        float simTime = exposure/timeDivision;
        float simDuration = nFrames*exposure; // in seconds
        int traceLength = (int) (simDuration/simTime);
        double sigmaPSF = 0.21*lambda/(NA*simPixelSize); // in simulated pixelsize
        float photonsPerFrame = photonFlux*exposure;

        IJ.log("kON: "+kON+" s^-1 (Tau_OFF: "+1000*(1/kON)+"ms)");
        IJ.log("kOFF: "+kOFF+" s^-1 (Tau_ON: "+1000*(1/kOFF)+"ms)");
        IJ.log("Photon flux: "+photonFlux+" photons/s");
        IJ.log("Acquisition time: "+simDuration+"s ("+nFrames+" frames at "+exposure*1000+"ms exposure)");

        IJ.log("PSF sigma: "+(sigmaPSF*simPixelSize)+" nm");

        boolean[][] moleculeTimeTraces = new boolean[nMolecules][traceLength];

        IJ.log("Simulating time trace from single molecule...");
        for (int i = 0; i < nMolecules; i++) {
            IJ.showProgress(i, nMolecules);
            boolean[] thisTrace = getTimeTraceSingleMolecule(kON, kOFF, simTime, simDuration);
            for (int j = 0; j < traceLength; j++) {
                moleculeTimeTraces[i][j] = thisTrace[j];
            }
        }

        int wSim = imp.getWidth();
        int hSim = imp.getHeight();
        int shrinkFactor = (int) (pixelSize/simPixelSize);

        ImageStack imsSim = new ImageStack(wSim, hSim);
        ImageStack imsFullSim = new ImageStack();
        ImageStack imsBg = new ImageStack();
        ImageStack imsSimNoBg = new ImageStack();
        float[] intensities;
        FloatProcessor fpSim;
        FloatProcessor fpSimBinned, fpBg;
        ImageProcessor ipBinned;
        float[] pixels;

        IJ.log("Generating fluorescence images...");

        for (int f = 0; f < nFrames; f++) {
            IJ.showProgress(f, nFrames);
            pixels = new float[wSim*hSim];

            intensities = new float[nMolecules];
            for (int m = 0; m < nMolecules; m++) {
                for (int t = 0; t < timeDivision; t++) {
                    if (moleculeTimeTraces[m][t + f*timeDivision]) intensities[m] += 1.0f/timeDivision;
                }
                pixels[xy[m][0]*wSim + xy[m][1]] += intensities[m]*photonsPerFrame; // ADC
            }

            fpSim = new FloatProcessor(wSim, hSim, pixels);
            fpSim.blurGaussian(sigmaPSF);
            if (showSimStack) imsSim.addSlice(fpSim);

            ipBinned = fpSim.duplicate();
            Binner binner = new Binner();
            ipBinned = binner.shrink(ipBinned, shrinkFactor, shrinkFactor, Binner.SUM);

            fpSimBinned = ipBinned.duplicate().convertToFloatProcessor();
//            fpSimBinned = fpSim.resize(w, h).convertToFloatProcessor();  // this doesn't do what we want. Doesn't add but decimates
            if (!disablePoissonNoise) fpSimBinned = addPoissonNoise(fpSimBinned).convertToFloatProcessor();
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

        ImagePlus impSimBinned = new ImagePlus("Simulated stack with binning", imsFullSim);
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


    public static boolean[] getTimeTraceSingleMolecule(float kON, float kOFF, float simTime, float duration){

        boolean showPlot = false;
        Random rnd = new Random();
        boolean ON;
        double p = rnd.nextDouble();
//        if (kON/kOFF > 1) ON = (p > (kOFF/kON)); // TODO: check that this is true
//        else ON = (p < (kON/kOFF));
        ON = (p < kON/(kON + kOFF));

//        IJ.log("START: Molecule is "+ON);

        int traceLength = (int) (duration/simTime);
        boolean[] switchingTrace = new boolean[traceLength];
        double[] switchingTraceDouble = null;
        if (showPlot) {
            switchingTraceDouble = new double[traceLength];
        }

//        IJ.log("Trace length: "+traceLength);

        int nextTimeChange;
        int traceID = 0;
        int nSwitches = 0;

        while ((traceID < traceLength)) {
//            IJ.log("Current trace ID: "+traceID);
//            IJ.log("Molecule is "+ON);
            p = rnd.nextDouble();
            if (ON) nextTimeChange = (int) (-Math.log(p) / (kOFF * simTime)); // in unit of simTime
            else nextTimeChange = (int) (-Math.log(p) / (kON * simTime)); // in unit of simTime
//            IJ.log("Next time change: "+(nextTimeChange*simTime)+" s");

            int i = 0;
            while ((i < nextTimeChange) && (traceID < traceLength)) {
                switchingTrace[traceID] = ON;
                if (showPlot) {
                    if (ON) switchingTraceDouble[traceID] = 1;
                    else switchingTraceDouble[traceID] = 0;
                }
                i++;
                traceID++;
            }
            ON = !ON;
            nSwitches ++;
        }


        if (showPlot) {
            Plot tracePlot = new Plot("Trace plot", "Time", "Switching state");
            tracePlot.add("line", switchingTraceDouble);
            tracePlot.show();

            // Calculating the ON and OFF times from the trace
            float avONtime = (float) sum(switchingTraceDouble)*simTime*1000/(nSwitches/2); // in ms
            float avOFFtime = 1000*(duration - (float) sum(switchingTraceDouble)*simTime)/(nSwitches/2);  // in ms
            IJ.log("Number of blinks: "+(nSwitches/2));
            IJ.log("Average ON time: "+avONtime+" ms");
            IJ.log("Average OFF time: "+avOFFtime+" ms");

        }

        //        IJ.log("Looping test");
//        Random rnd = new Random();
//        double p;
//        for (int i = 0; i < 10; i++) {
//            p = rnd.nextDouble();
//            double nextTimeChange = (-Math.log(p) / kOFF); // in unit of simTime
//            IJ.log("p="+p);
//            IJ.log("Next time: "+1000*nextTimeChange+" ms");
//        }

        return switchingTrace;
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
