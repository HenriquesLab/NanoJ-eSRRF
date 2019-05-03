package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import nanoj.liveSRRF.CSVImport;
import nanoj.liveSRRF.NNDistance_CL;
import org.scijava.table.DefaultFloatTable;

public class NNAnalysis_ implements PlugIn {

    private int nLocs;

    // TODO: how to ignore the dark areas? BG localizations? run DBSCAN on localization on TS beforehand!


    private DefaultFloatTable dataTable;

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        OpenDialog dialog = new OpenDialog("Please choose a localization file"); // TODO: handle cancel
        String filePath = dialog.getPath();
        try {
            IJ.log("Reading file...");
            CSVImport csvImport = new CSVImport(filePath);
            dataTable = csvImport.dataTable;
        } catch (Exception e) {
            IJ.log("Computer says no.");
        }

        nLocs = dataTable.getRowCount();
        int nFrames = Math.round(dataTable.get(0, nLocs - 1));
        IJ.log("Number of localizations: " + nLocs);
        IJ.log("Number of columns: " + dataTable.getColumnCount());
        IJ.log("Number of frames: " + nFrames);


        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("NanoJ - Localisation density analysis");
        gd.addNumericField("Size of temporal blocks (frames): ", 100, 0);
        gd.addNumericField("Radius for density calculation (nm): ", 100, 2);
        gd.addNumericField("Number of spatial blocks per axis: ", 4, 0);
        gd.addNumericField("Image width/height (nm): ", 40960, 2);
        gd.addCheckbox("Show histograms ", false);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int blockSize = (int) gd.getNextNumber();
        float densityRadius = (float) gd.getNextNumber();
        int blockPerAxis = (int) gd.getNextNumber();
        float imageSize = (float) gd.getNextNumber();
        boolean showHistogram = gd.getNextBoolean();


        NNDistance_CL nndCalculator = new NNDistance_CL(densityRadius, blockPerAxis, imageSize);

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

            int nStart = getFramePosition(ntB * blockSize + 1) + 1;
            int nEnd = getFramePosition((ntB + 1) * blockSize + 1);
            nLocsPerBlock[ntB] = nEnd - nStart + 1;

            IJ.log("Block #" + (ntB + 1) + "/" + nTempBlocks + " (using localizations #" + nStart + " to #" + nEnd + ")");

            if (nLocsPerBlock[ntB] > 1) {
                float[] xyDataArray = new float[2 * nLocsPerBlock[ntB]];
                for (int i = 0; i < nLocsPerBlock[ntB]; i++) {
                    xyDataArray[i] = dataTable.get(1, i + nStart);
                    xyDataArray[i + nLocsPerBlock[ntB]] = dataTable.get(2, i + nStart);
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
        double[] globalMeanStd = getMeanStdFromArray(meanNNDarray);
        IJ.log("Temporal block size: " + blockSize + " frames.");
        IJ.log(String.format("NN distance (nm): %.2f +/- %.2f", globalMeanStd[0], globalMeanStd[1]));

        IJ.log("-----------------------------");
        IJ.log("All done. Cheerio !");

    }


    private int getFramePosition(int frameNumber){

        int locPosition = -1;
        for (int i=0; i<nLocs; i++){
            if (dataTable.get(0, i) < frameNumber) locPosition = i;
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
