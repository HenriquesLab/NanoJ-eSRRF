package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.CSVImport;
import nanoj.liveSRRF.NNDistance_CL;
import org.scijava.table.DefaultFloatTable;

public class NNAnalysis_ implements PlugIn {

    private int nLocs;
    private int nColX, nColY, nColFrame, nFrameToChop;
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());



    // TODO: how to ignore the dark areas? BG localizations? run DBSCAN on localization on TS beforehand!


    private DefaultFloatTable dataTable;

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");


        OpenDialog dialog = new OpenDialog("Please choose a localization file"); // TODO: handle cancel
        String filePath = dialog.getPath();

        String[] simChoice = new String[2];
        simChoice[0] = "From ThunderSTORM";
        simChoice[1] = "From SureSim";
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("NanoJ - Localisation density analysis");
        gd.addChoice("Loc file type: ", simChoice, prefs.get("simType", simChoice[0]));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        String simType = gd.getNextChoice();

        if (simType.equals(simChoice[0])) {
            try {
                IJ.log("Reading file...");
                IJ.log(filePath);
                CSVImport csvImport = new CSVImport(filePath, ",");
                dataTable = csvImport.dataTable;
            } catch (Exception e) {
                IJ.log("Computer says no.");
            }
            nColX = 1;
            nColY = 2;
            nColFrame = 0;
        }

        else{
            try {
                IJ.log("Reading file...");
                IJ.log(filePath);
                CSVImport csvImport = new CSVImport(filePath, " ");
                dataTable = csvImport.dataTable;
            } catch (Exception e) {
                IJ.log("Computer says no.");
            }
            nColX = 0;
            nColY = 1;
            nColFrame = 3;
        }

        nLocs = dataTable.getRowCount();
        int nFrames = Math.round(dataTable.get(nColFrame, nLocs - 1));
        IJ.log("Number of localizations: " + nLocs);
        IJ.log("Number of columns: " + dataTable.getColumnCount());
        IJ.log("Number of frames: " + nFrames);

        gd = new NonBlockingGenericDialog("NanoJ - Localisation density analysis");
        gd.addNumericField("Size of temporal blocks (frames): ", prefs.get("blockSize", 100), 0);
        gd.addNumericField("Radius for density calculation (nm): ", prefs.get("densityRadius", 100), 2);
        gd.addNumericField("Number of spatial blocks per axis: ", prefs.get("blockPerAxis", 4), 0);
        gd.addNumericField("Image width/height (nm): ", prefs.get("imageSize", 40960), 2);
        gd.addNumericField("Number of frames to ignore (at start & end): ", prefs.get("nFrameToChop", 10), 0);
        gd.addNumericField("GPU block size: ", prefs.get("blockLength", 20000), 0);
        gd.addCheckbox("Show histograms ", prefs.get("showHistogram", false));
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int blockSize = (int) gd.getNextNumber();
        float densityRadius = (float) gd.getNextNumber();
        int blockPerAxis = (int) gd.getNextNumber();
        float imageSize = (float) gd.getNextNumber();
        nFrameToChop = (int) gd.getNextNumber();
        int blockLength = (int) gd.getNextNumber();
        boolean showHistogram = gd.getNextBoolean();


        prefs.set("simType", simType);
        prefs.set("blockSize", blockSize);
        prefs.set("densityRadius", densityRadius);
        prefs.set("blockPerAxis", blockPerAxis);
        prefs.set("imageSize", imageSize);
        prefs.set("nFrameToChop", nFrameToChop);
        prefs.set("blockLength", blockLength);
        prefs.set("showHistogram", showHistogram);


        nFrames -= 2*nFrameToChop; // adjust the number of frames

        NNDistance_CL nndCalculator = new NNDistance_CL(densityRadius, blockPerAxis, imageSize, blockLength);

        int nTempBlocks = nFrames / blockSize;
        IJ.log("Number of blocks to compute: " + nTempBlocks);
        IJ.log("Radius for density calculation (nm): " + densityRadius);
        IJ.log("Number of spatial blocks per axis: " + blockPerAxis);

        int[] nLocsPerBlock = new int[nTempBlocks];
        double[] meanNNDarray = new double[nTempBlocks];
//        double[] stdNNDarray = new double[nTempBlocks];
        double[] tempMeanStdNND;

        double[] meanDensityArray = new double[nTempBlocks];
//        double[] stdDensityArray = new double[nTempBlocks];
        double[] tempMeanStdDensity;

        for (int ntB = 0; ntB < nTempBlocks; ntB++) {

            int nStart = getFramePosition(ntB * blockSize + 1 + nFrameToChop) + 1;
            int nEnd = getFramePosition((ntB + 1) * blockSize + 1 + nFrameToChop);
            nLocsPerBlock[ntB] = nEnd - nStart + 1;

            IJ.log("Block #" + (ntB + 1) + "/" + nTempBlocks + " (using localizations #" + nStart + " to #" + nEnd + ")");

            if (nLocsPerBlock[ntB] > 1) {
                float[] xyDataArray = new float[2 * nLocsPerBlock[ntB]];
                for (int i = 0; i < nLocsPerBlock[ntB]; i++) {
                    xyDataArray[i] = dataTable.get(nColX, i + nStart);
                    xyDataArray[i + nLocsPerBlock[ntB]] = dataTable.get(nColY, i + nStart);
                }

                nndCalculator.calculateNND(xyDataArray);
                double[] buffersArray = nndCalculator.readBuffers();
                nndCalculator.release();

                double[] nndArray = new double[nLocsPerBlock[ntB]];
                double[] densityArray = new double[nLocsPerBlock[ntB]];

                for (int i=0; i<nLocsPerBlock[ntB]; i++){
                    nndArray[i] = buffersArray[i];
                    densityArray[i] = buffersArray[i+nLocsPerBlock[ntB]];
                }

                if (showHistogram) {
                    // Plot the histogram
                    IJ.log("Plotting values...");
                    Plot valuePlot = new Plot("NN distance values", "Localization #", "NN distance (nm)");
                    valuePlot.add("circle", nndArray);
                    valuePlot.show();

                    // Plot the histogram
                    IJ.log("Plotting histogram...");
                    Plot histPlot = new Plot("NN distance histogram", "NN distance (nm)", "Occurences");
                    histPlot.addHistogram(nndArray);
                    histPlot.show();
                }

                tempMeanStdNND = getMeanStdFromArray(nndArray);
                tempMeanStdDensity = getMeanStdFromArray(densityArray);

                meanNNDarray[ntB] = tempMeanStdNND[0];
//                stdNNDarray[ntB] = tempMeanStdNND[1];
                meanDensityArray[ntB] = tempMeanStdDensity[0];
//                stdDensityArray[ntB] = tempMeanStdDensity[1];

            } else {
                IJ.log("No localizations found in these frames...");
            }

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }

            IJ.showProgress((double) (ntB + 1) / (double) nTempBlocks);

        }


        Plot nndBlockPlot = new Plot("NN distance vs. temporal block", "Temporal block #", "Mean NND (nm)");
        nndBlockPlot.add("circle", meanNNDarray);
        nndBlockPlot.show();

        double[] nLocsPerBlockDouble = new double[nTempBlocks];
        for (int i = 0; i < nTempBlocks; i++) {
            nLocsPerBlockDouble[i] = (double) nLocsPerBlock[i];
        }

        Plot nLocsBlockPlot = new Plot("Number of locs vs. temporal block", "Temporal block #", "Number of locs per block");
        nLocsBlockPlot.add("circle", nLocsPerBlockDouble);
        nLocsBlockPlot.show();

        Plot densityBlockPlot = new Plot("Density vs. temporal block", "Temporal block #", "Density (#loc/defined area)");
        densityBlockPlot.add("circle", meanDensityArray);
        densityBlockPlot.show();

        ImagePlus impNND = new ImagePlus("NNd map", nndCalculator.imsNND);
        impNND.show();
        IJ.run(impNND, "SQUIRREL-FRC", "");
        ImagePlus impDensity = new ImagePlus("Density map", nndCalculator.imsDensity);
        impDensity.show();
        IJ.run(impDensity, "SQUIRREL-FRC", "");


        IJ.log("-----------------------------");
        double[] globalMeanStdNND = getMeanStdFromArray(meanNNDarray);
        IJ.log("Temporal block size: " + blockSize + " frames.");
        IJ.log(String.format("NN distance (nm): %.2f +/- %.2f", globalMeanStdNND[0], globalMeanStdNND[1]));

        double[] globalMeanStdDensity = getMeanStdFromArray(meanDensityArray);
        double area = Math.PI*densityRadius*densityRadius/(1e6);
        globalMeanStdDensity[0] /= area;
        globalMeanStdDensity[1] /= area;

        IJ.log(String.format("Density (#loc/um^2): %.2f +/- %.2f (using %.2f nm as a radius)", globalMeanStdDensity[0], globalMeanStdDensity[1], densityRadius));

        IJ.log("-----------------------------");
        IJ.log("All done. Cheerio !");

    }


    private int getFramePosition(int frameNumber){

        int locPosition = -1;
        for (int i=0; i<nLocs; i++){
            if (dataTable.get(nColFrame, i) < frameNumber) locPosition = i;
        }

        return locPosition;
    }

    private double[] getMeanStdFromArray(double[] array){

        int arraySize = array.length;
        double[] meanStd = new double[2];
        for (double v : array) {
            meanStd[0] += v;
            meanStd[1] += v * v;
        }

        meanStd[0] /= arraySize;
        meanStd[1] /= arraySize;
        meanStd[1] = Math.sqrt(meanStd[1] - meanStd[0] * meanStd[0]);

        return meanStd;
    }

}
