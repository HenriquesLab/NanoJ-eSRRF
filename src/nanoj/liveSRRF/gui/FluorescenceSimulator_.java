package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.util.LinkedHashMap;
import java.util.Map;

import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import java.util.Random;



public class FluorescenceSimulator_ implements PlugIn {

    ImagePlus imp;

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

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Fluorescence simulator");
        gd.addNumericField("Emission wavelength (nm):", 700, 2);
        gd.addNumericField("NA:", 1.15, 2);
        gd.addNumericField("kON (s^-1):", 10, 2);
        gd.addNumericField("kOFF (s^-1):", 10, 2);
        gd.addNumericField("Number of frames: ", 100, 0);
        gd.addNumericField("Exposure time (ms): ", 10, 1);
        gd.addNumericField("Pixel size (nm): ", 100, 1);
        gd.showDialog();

        float lambda = (float) gd.getNextNumber();
        float NA = (float) gd.getNextNumber();
        float kON = (float) gd.getNextNumber();
        float kOFF = (float) gd.getNextNumber();
        float exposure = ((float) gd.getNextNumber())/1000; // convert to seconds
        int nFrames = (int) gd.getNextNumber();
        float pixelSize = (float) gd.getNextNumber();

        Calibration cal = imp.getCalibration();
        if (cal.pixelHeight != cal.pixelWidth) return;
        float simPixelSize = (float) cal.pixelHeight;

        ImageStack imsGT = imp.getStack();
        int[][] xy = getMoleculePositionsFromGroundTruth(imsGT);
        int nMolecules = xy.length;

        int timeDivision = 10;
        float simTime = exposure/timeDivision;
        float simDuration = nFrames*exposure;
        int traceLength = (int) (simDuration/simTime);
        float sigmaPSF = 0.21f*lambda/(NA*simPixelSize); // in simulated pixelsize

        IJ.log("PSF sigma: "+(sigmaPSF*simPixelSize)+" nm");

        boolean[][] moleculeTimeTraces = new boolean[nMolecules][traceLength];

        IJ.log("Simulating time trace from single molecule...");

        for (int i = 0; i < nMolecules; i++) {
            boolean[] thisTrace = getTimeTraceSingleMolecule(kON, kOFF, simTime, simDuration);
            for (int j = 0; j < traceLength; j++) {
                moleculeTimeTraces[i][j] = thisTrace[j];
            }
        }

        int wSim = imp.getWidth();
        int hSim = imp.getHeight();
        int w = (int) (wSim*pixelSize/simPixelSize);
        int h = (int) (hSim*pixelSize/simPixelSize);

        ImageStack imsSim = new ImageStack(wSim, hSim);
        float[] intensities;
        FloatProcessor fpSim;

        for (int f = 0; f < nFrames; f++) {
            fpSim = new FloatProcessor(wSim, hSim, new float[wSim*hSim]);
            intensities = new float[nMolecules];
            for (int m = 0; m < nMolecules; m++) {
                for (int t = 0; t < timeDivision; t++) {
                    if (moleculeTimeTraces[m][t + f*timeDivision]) intensities[m] += 1.0f/timeDivision;
                }
                fpSim.setf(xy[m][0], xy[m][1], intensities[m]);
            }
            imsSim.addSlice(fpSim);
        }

        ImagePlus impSim = new ImagePlus("Simulated stack", imsSim);
        impSim.show();

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
//        IJ.log("Number of molecules: "+nMolecules);
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
        if (kON/kOFF > 1) ON = (p > (kOFF/kON));
        else ON = (p < (kON/kOFF));
//        IJ.log("START: Molecule is "+ON);

        int traceLength = (int) (duration/simTime);

        boolean[] switchingTrace = new boolean[traceLength];
        double[] switchingTraceDouble = null;
        if (showPlot) {
            switchingTraceDouble = new double[traceLength];
        }

//        IJ.log("Trace length: "+traceLength);

        switchingTrace[0] = ON;

        int nextTimeChange;
        int traceID = 0;

//        int exitDoorID = 0; // there to ensure that we exit the while loop
//        while ((traceID < traceLength) && (exitDoorID<10)) {

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
//            exitDoorID++;
        }

        if (showPlot) {
            Plot tracePlot = new Plot("Trace plot", "Time", "Switching state");
            tracePlot.add("line", switchingTraceDouble);
            tracePlot.show();
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

}
