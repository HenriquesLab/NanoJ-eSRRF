package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import nanoj.liveSRRF.CSVImport;
import nanoj.liveSRRF.NNDistance_CL;
import org.scijava.table.DefaultFloatTable;

public class NNAnalysis_ implements PlugIn {

    int nLocs;

    // TODO: - correct for the edges? ignore localizations that are too close to the edges?
    // TODO: compute in spatial blocks
    // TODO: how to ignore the dark areas? BG localizations? run DBSCAN on localization on TS beforehand!


    public DefaultFloatTable dataTable;

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        OpenDialog dialog = new OpenDialog("Please choose a localization file"); // TODO: handle cancel
        String filePath = dialog.getPath();
        try{
            IJ.log("Reading file...");
            CSVImport csvImport = new CSVImport(filePath);
            dataTable = csvImport.dataTable;}
        catch(Exception e){
            IJ.log("Computer says no.");
        }

        nLocs = dataTable.getRowCount();
        int nFrames = Math.round(dataTable.get(0,nLocs-1));
        IJ.log("Number of localizations: "+nLocs);
        IJ.log("Number of columns: "+dataTable.getColumnCount());
        IJ.log("Number of frames: "+nFrames);


        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("NanoJ - Nearest-neighbour analysis");
        gd.addNumericField("Size of temporal blocks (frames): ", 100, 0);
        gd.addNumericField("Radius for density calculation (nm): ", 100, 2);
        gd.addCheckbox("Show histograms ", false);
        gd.showDialog();
        if (gd.wasCanceled()) return;

        int blockSize = (int) gd.getNextNumber();
        float densityRadius = (float) gd.getNextNumber();
        boolean showHistogram = gd.getNextBoolean();


        NNDistance_CL nndCalculator = new NNDistance_CL(densityRadius);
        int nBlocks = nFrames/blockSize;
        IJ.log("Number of blocks to compute: "+nBlocks);
        IJ.log("Radius for density calculation (nm): "+densityRadius);

        int[] nLocsPerBlock = new int[nBlocks];
        double[] meanNNDarray = new double[nBlocks];
        double[] stdNNDarray = new double[nBlocks];
        double[] tempMeanStdNND;

        double[] meanDensityArray = new double[nBlocks];
        double[] stdDensityArray = new double[nBlocks];
        double[] tempMeanStdDensity;

        for (int nB = 0; nB < nBlocks; nB++) {

            int nStart = getFramePosition(nB*blockSize+1)+1;
            int nEnd = getFramePosition((nB+1)*blockSize+1);
            nLocsPerBlock[nB] = nEnd-nStart+1;

            IJ.log("Block #" + (nB + 1) + "/" + nBlocks + " (using localizations #" + nStart + " to #" + nEnd + ")");

            if (nLocsPerBlock[nB]>1) {

                float[] xyDataArray = new float[2 * nLocsPerBlock[nB]];
                for (int i = 0; i < nLocsPerBlock[nB]; i++) {
                    xyDataArray[i] = dataTable.get(1, i + nStart);
                    xyDataArray[i + nLocsPerBlock[nB]] = dataTable.get(2, i + nStart);
                }

                nndCalculator.calculateNND(xyDataArray);
                double[] nndArray = nndCalculator.readNNDbuffer();
                int[] densityArray = nndCalculator.readDensitybuffer();
                nndCalculator.release();

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

                meanNNDarray[nB] = tempMeanStdNND[0];
                stdNNDarray[nB] = tempMeanStdNND[1];
                meanDensityArray[nB] = tempMeanStdDensity[0];
                stdDensityArray[nB] = tempMeanStdDensity[1];
            }
            else{
                IJ.log("No localizations found in these frames...");
            }

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }

            IJ.showProgress((double) (nB+1) / (double) nBlocks);

        }


        Plot nndBlockPlot = new Plot("NN distance vs. temporal block", "Temporal block #", "Mean NND (nm)");
        nndBlockPlot.add("circle", meanNNDarray);
        nndBlockPlot.show();

        double[] nLocsPerBlockDouble = new double[nBlocks];
        for (int i=0; i<nBlocks; i++){
            nLocsPerBlockDouble[i] = (double) nLocsPerBlock[i];
        }

        Plot nLocsBlockPlot = new Plot("Number of locs vs. temporal block", "Temporal block #", "Number of locs per block");
        nLocsBlockPlot.add("circle", nLocsPerBlockDouble);
        nLocsBlockPlot.show();

        Plot densityBlockPlot = new Plot("Density vs. temporal block", "Temporal block #", "Density (#loc/defined area)");
        densityBlockPlot.add("circle", meanDensityArray);
        densityBlockPlot.show();

        IJ.log("-----------------------------");
        double[] globalMeanStd = getMeanStdFromArray(meanNNDarray);
        IJ.log("Temporal block size: "+blockSize+" frames.");
        IJ.log(String.format("NN distance (nm): %.2f +/- %.2f", globalMeanStd[0], globalMeanStd[1]));

        IJ.log("-----------------------------");
        IJ.log("All done. Cheerio !");

    }


    public int getFramePosition(int frameNumber){

        int locPosition = -1;
        for (int i=0; i<nLocs; i++){
            if (dataTable.get(0, i) < frameNumber) locPosition = i;
        }

        return locPosition;
    }

    public double[] getMeanStdFromArray(double[] array){

        int arraySize = array.length;
        double[] meanStd = new double[2];
        for (int i = 0; i < arraySize; i++) {
            meanStd[0] += array[i];
            meanStd[1] += array[i] * array[i];
        }

        meanStd[0] /= arraySize;
        meanStd[1] /= arraySize;
        meanStd[1] = Math.sqrt(meanStd[1] - meanStd[0] * meanStd[0]);

        return meanStd;
    }

    public double[] getMeanStdFromArray(int[] array){

        int arraySize = array.length;
        double[] meanStd = new double[2];
        for (int i = 0; i < arraySize; i++) {
            meanStd[0] += (double) array[i];
            meanStd[1] += (double) array[i] * array[i];
        }

        meanStd[0] /= arraySize;
        meanStd[1] /= arraySize;
        meanStd[1] = Math.sqrt(meanStd[1] - meanStd[0] * meanStd[0]);

        return meanStd;
    }


}
