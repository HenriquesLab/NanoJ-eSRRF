package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.measure.ResultsTable;
import ij.plugin.Concatenator;
import ij.plugin.PlugIn;
import nanoj.core.java.io.LoadNanoJTable;

import java.io.IOException;
import java.util.Map;

import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.liveSRRF.GetShiftAndTiltRCCM.applyCorrection;
import static nanoj.liveSRRF.gui.GetSpatialCalibrationMFMdata_.getMaxfromArray;

public class ApplyCalibrationMFMdata_ implements PlugIn {

    private double[] shiftX, shiftY, theta;
    private final String MFMApplyCalibVersion = "v0.1";

    @Override
    public void run(String s) {
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("MFM Calibration ("+MFMApplyCalibVersion+")");

        ImageStack ims = imp.getImageStack();
        int width = imp.getWidth();
        int height = imp.getHeight();
        int nFrames = imp.getStackSize();

        IJ.log("Getting calibration file...");
        String calibTablePath = IJ.getFilePath("Choose Drift-Table to load...");
        Map<String, double[]> calibTable;
        try {
            calibTable = new LoadNanoJTable(calibTablePath).getData();
            shiftX = calibTable.get("X-shift (pixels)");
            shiftY = calibTable.get("Y-shift (pixels)");
            theta = calibTable.get("Theta (degrees)");
            ResultsTable rt = dataMapToResultsTable(calibTable);
            rt.show("Calibration-Table");
        } catch (IOException e) {
            IJ.log("Catching exception...");
            e.printStackTrace();
        }

        int nSplits = (int) Math.sqrt(shiftX.length);
        int cropSizeX = Math.round(width/nSplits); // TODO: assume that size in X and Y are the same?
        int cropSizeY = Math.round(height/nSplits); // TODO: assume that size in X and Y are the same?

        ImagePlus[] impCorrected = new ImagePlus[nSplits*nSplits];

        int x,y;
        double[] shiftXslice = new double[nFrames];
        double[] shiftYslice = new double[nFrames];
        double[] thetaSlice = new double[nFrames];

        int[] offsetMFMarray;
//        if (nSplits == 3){
        offsetMFMarray = new int[] {2,3,4,-1,0,1,-4,-3,-2};
        // TODO: add other cases where it's not 9 positions
//        }
        int maxOffset = getMaxfromArray(offsetMFMarray);

        IJ.log("Applying calibration...");
        for (int j = 0; j < nSplits; j++){
            for (int i = 0; i < nSplits; i++){
                for (int k = 0; k < nFrames; k++){
                    shiftXslice[k] = shiftX[j*nSplits+i];
                    shiftYslice[k] = shiftY[j*nSplits+i];
                    thetaSlice[k] = theta[j*nSplits+i];
                }
                IJ.log("X-shift: "+shiftXslice[0]);

                x = Math.round((width/nSplits)*i);
                y = Math.round((height/nSplits)*j);
                ImageStack imsTemp = ims.crop(x, y, 0, cropSizeX,cropSizeY, nFrames);
                ImagePlus impTemp = new ImagePlus("Substack "+i+"/"+j,imsTemp);
                impCorrected[offsetMFMarray[j*nSplits+i]+maxOffset] = applyCorrection(impTemp, shiftXslice, shiftYslice, thetaSlice);
            }
        }

        IJ.log("Reshaping data...");
        Concatenator concatenator = new Concatenator(); // TODO: convert this into the right stack dimensions
        ImagePlus impStack = concatenator.concatenate(impCorrected, false);
        impStack.show();

        IJ.log("------------");
        IJ.log("All done.");
    }
}
