package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.gui.Plot;
import ij.io.OpenDialog;
import ij.plugin.PlugIn;
import nanoj.liveSRRF.CSVImport;
import nanoj.liveSRRF.NNDistance_CL;
import org.scijava.table.DefaultFloatTable;

public class NNAnalysis_ implements PlugIn {

    public DefaultFloatTable dataTable;

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        OpenDialog dialog = new OpenDialog("Please choose a localization file");
        String filePath = dialog.getPath();
        try{
            IJ.log("Reading file...");
            CSVImport csvImport = new CSVImport(filePath);
            dataTable = csvImport.dataTable;}
        catch(Exception e){
            IJ.log("Computer says no.");
        }

        int nLocs = dataTable.getRowCount();
        IJ.log("Number of localizations: "+nLocs);
        IJ.log("Number of columns: "+dataTable.getColumnCount());

        float[] xyDataArray = new float[2*nLocs];
        for (int i = 0; i < nLocs; i++){
            xyDataArray[i] = dataTable.get(1,i);
            xyDataArray[i+nLocs] = dataTable.get(2,i);
        }

        NNDistance_CL nndCalculator = new NNDistance_CL(xyDataArray);
        double[] nndArray = nndCalculator.readNNDbuffer();

        nndCalculator.release();

        // Plot the histogram
        IJ.log("Plotting values...");
        Plot valuePlot = new Plot("NN distance values","NN distance (nm)", "Occurences");
        valuePlot.add("circle", nndArray);
        valuePlot.show();

        // Plot the histogram
        IJ.log("Plotting histogram...");
        Plot histPlot = new Plot("NN distance histogram","NN distance (nm)", "Occurences");
        histPlot.addHistogram(nndArray);
        histPlot.show();

        IJ.log("------------- Statistics -----------");
        double meanNND = 0;
        double stdNND = 0;
        for (int i=0; i<nLocs; i++){
            meanNND += nndArray[i]/nLocs;
            stdNND += nndArray[i]*nndArray[i]/nLocs;
        }
        stdNND = Math.sqrt(stdNND - meanNND*meanNND);

        IJ.log(String.format("Mean: %.2f nm", meanNND));
        IJ.log(String.format("STD: %.2f nm", stdNND));

        IJ.log("-----------------------------");
        IJ.log("All done. Cheerio !");

    }
}
