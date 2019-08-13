package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core.java.featureExtraction.Peaks;
import nanoj.core2.NanoJPrefs;
import nanoj.core2.NanoJProfiler;
import nanoj.liveSRRF.Fit1DGaussian;
import org.python.modules.math;

import java.awt.*;


public class GetDefocus_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private NanoJProfiler prof = new NanoJProfiler();
    int width, height, nSlices;


    public void run(String s) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000*imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep*100d)/100d;

        // initilizaing string for focus estimation choice
//        String[] FocusMethod = new String[2];
//        FocusMethod[0] = "Intensity max";
//        FocusMethod[1] = "Entropy";

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Get defocus from Image splitter:");
        gd.addNumericField("Number of split/axis:", prefs.get("nSplits", 3), 0);
//        gd.addChoice("Focus position estimation method", FocusMethod, FocusMethod[0]);

        gd.addNumericField("Maximum number of beads to analyse:", prefs.get("nPeaks", 5000), 0);
        gd.addNumericField("Crop level:", prefs.get("cropLevel", 0.5f), 2);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int nSplits = (int) gd.getNextNumber();
        int nPeaks = (int) gd.getNextNumber();
        float cropLevel = (float) gd.getNextNumber();

        prefs.set("nSplits", (float) nSplits);
        prefs.set("nPeaks", (float) nPeaks);
        prefs.set("cropLevel", (float) cropLevel);

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("z-step: "+zStep+" nm");
        IJ.log("Level for cropping: "+cropLevel);

        ImageStack ims = imp.getImageStack();
        this.nSlices = ims.getSize();
        this.width = ims.getWidth();
        this.height = ims.getHeight();

        Peaks pks = new Peaks();
        int radius = 3; // pixels TODO: make this adjustable for larger blobs?
        float maxOverlap = 0.0f; // TODO: make this adjustable?

        double[] zPos = new double[nSlices];
        for (int n = 0; n < nSlices; n++) zPos[n] = n*zStep;

        Plot defocusPlot = new Plot("Defocus plot", "Z (nm)", "Intensity (AU)");
//        Plot EntropyPlot = new Plot("Entropy plot", "Z (nm)", "Entropy (AU)");

        Fit1DGaussian fitting;
        float[] defocusArray;
        float[][] arrayPeaks;
        float[] fitResults;
        double[][] modelArray;
        double[] focusPositionArray = new double[nSplits];
//        ShannonEntropyMeasure em;
//        double[] EntropyArray;



        for (int n = 0; n < nSplits; n++){

            IJ.log("------------ Slice "+(n+1)+"/"+nSplits+" ----------");
            ImageStack imsSlice = ims.crop(0, (int) math.floor(n*height / nSplits), 0, width, (int) math.floor(height / nSplits), nSlices).convertToFloat();
//            IJ.log("Size of the cropped images: "+imsSlice.getWidth()+" x "+imsSlice.getHeight());

            int id = prof.startTimer();
            ImageStack imsSum = sumImageStack(imsSlice);
            prof.recordTime("Stack summation: ", prof.endTimer(id));

            id = prof.startTimer();
            arrayPeaks = pks.getPeaks(imsSum.getProcessor(1), nPeaks, radius, maxOverlap);
            IJ.log("Number of peaks found: " + arrayPeaks.length);
            prof.recordTime("Get peaks: ", prof.endTimer(id));

            id = prof.startTimer();
            defocusArray = getDefocusArray(arrayPeaks, imsSlice);
            fitting = new Fit1DGaussian(defocusArray);
            fitting.cropDataArray(cropLevel);
            fitResults = fitting.calculate();
            focusPositionArray[n] = fitResults[1]*zStep;
            prof.recordTime("Defocus array: ", prof.endTimer(id));


            IJ.log("--- Initial guesses results ---");
            IJ.log("Initial x0: "+(fitting.initX0)*zStep+" nm");
            IJ.log("Initial Sigma: "+(fitting.initSigma)*zStep+" nm");

            IJ.log("--- Fit results ---");
            IJ.log("Amp: "+fitResults[0]);
            IJ.log("x0: "+(fitResults[1])*zStep+" nm");
            IJ.log("Sigma: "+(fitResults[2])*zStep+" nm");
            IJ.log("BG: "+fitResults[3]);

            modelArray = fitting.fittedCurve();
            defocusPlot.setColor(Color.black);
            defocusPlot.add("line", zPos, modelArray[0]);
            defocusPlot.setColor(Color.red);
            defocusPlot.add("line", zPos, modelArray[1]);

//            em = new ShannonEntropyMeasure(imsSlice);
//            EntropyArray = em.EntropyArray;
//
//            EntropyPlot.add("line", zPos, EntropyArray);


        }

        defocusPlot.setLimitsToFit(true);
        defocusPlot.show();
//        EntropyPlot.show();
        IJ.log("--- Profiler report ---");
        IJ.log(prof.report()); // TODO: look into profiler since iterations add up between consective runs (by design?)


        double dz;
        IJ.log("--- DEFOCUS RESULTS ---");
        for (int n = 0; n < nSplits; n++) {
            dz = focusPositionArray[n]-focusPositionArray[0];
            dz = (double) Math.round(dz*10d)/10d;
            IJ.log("Slice #"+n+": "+dz+" nm");
        }
    }


    // Support functions -----------------------------------------------------------------------------------------------
    public ImageStack sumImageStack(ImageStack ims) {

        float[] dataSum = new float[width*height];

        for (int n = 1; n < nSlices; n++) {
            float[] data = (float[]) ims.getPixels(n);
            for (int p=0; p<data.length; p++) dataSum[p] += data[p]/nSlices;
        }

        ImageStack imsSum = new ImageStack(width, height);
        imsSum.addSlice(new FloatProcessor(width, height, dataSum));
        return imsSum;
    }


    public float[] getDefocusArray(float[][] arrayPeaks, ImageStack ims) {

        int nPeaks = arrayPeaks.length;
        float[] defocusArray = new float[nSlices];

        for (int n = 0; n < nSlices; n++) {
            for (int p = 0; p < nPeaks; p++) {
//                int offset = (int) arrayPeaks[p][2] * width + (int) arrayPeaks[p][1]; // TODO: use getPixels
//                defocusArray[n] += ((float[]) ims.getPixels(n))[offset];
                defocusArray[n] += (float) (ims.getVoxel((int) arrayPeaks[p][1], (int) arrayPeaks[p][2], n))/nPeaks;
            }
        }

        return defocusArray;

    }
}