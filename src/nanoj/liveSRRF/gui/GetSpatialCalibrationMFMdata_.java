package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.io.FileInfo;
import ij.measure.ResultsTable;
import ij.plugin.Concatenator;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.GetShiftAndTiltRCCM;

import java.io.IOException;
import java.util.LinkedHashMap;
import java.util.Map;

import static nanoj.core.java.imagej.FilesAndFoldersTools.getSavePath;
import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.core.java.io.SaveNanoJTable.saveNanoJTable;
import static nanoj.liveSRRF.GetShiftAndTiltRCCM.applyCorrection;
import static nanoj.liveSRRF.StackProjections.calculateAverage;

public class GetSpatialCalibrationMFMdata_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private final String MFMGetCalibVersion = "v0.1";


    public void run(String s) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000 * imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep * 100d) / 100d;

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Spatial calibration from MFM data ("+MFMGetCalibVersion+")");
        gd.addNumericField("Number of split/axis:", prefs.get("nSplits", 3), 0);
        gd.addNumericField("Axial spacing (in nm):", prefs.get("axialSpacing", 390), 1);
        gd.addNumericField("Border crop (in pixels):", prefs.get("borderCrop", 10), 0);
        gd.addMessage("-=-= Search parameters =-=-\n");
        gd.addNumericField("Number of angles:", prefs.get("nAngles", 20), 0);
        gd.addNumericField("Max. angle (in degrees):", prefs.get("maxAngle", 2), 1);
        gd.addCheckbox("Apply correction", prefs.get("applyCorr", true));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        String header = "Choose where to save Calibration table...";
        String filePath = getSavePath(header, imp.getTitle(), ".njt");
        if (filePath == null) return;
        filePath = filePath.replace(".njt", "");

        int nSplits = (int) gd.getNextNumber();
        float axialSpacing = (float) gd.getNextNumber();
        int borderCrop = (int) gd.getNextNumber();

        int nAngles = (int) gd.getNextNumber();
        float maxAngle = (float) gd.getNextNumber();
        boolean applyCorr = gd.getNextBoolean();

        prefs.set("nSplits", (float) nSplits);
        prefs.set("axialSpacing", axialSpacing);
        prefs.set("borderCrop", (float) borderCrop);
        prefs.set("nAngles", (float) nAngles);
        prefs.set("maxAngle", maxAngle);
        prefs.set("applyCorr", applyCorr);


        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("MFM Calibration "+MFMGetCalibVersion);
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("Axial scan step size (in nm): "+zStep);

        int width, height, nSlices;
        ImageStack ims = imp.getImageStack();
        nSlices = ims.getSize();
        width = ims.getWidth();
        height = ims.getHeight();

        int[] offsetMFMarray;
//        if (nSplits == 3){
        offsetMFMarray = new int[] {2,3,4,-1,0,1,-4,-3,-2};
        // TODO: add other cases where it's not 9 positions
//        }

        int cropSizeX = Math.round(width/nSplits) - 2*borderCrop; // TODO: assume that size in X and Y are the same?
        int cropSizeY = Math.round(height/nSplits) - 2*borderCrop; // TODO: assume that size in X and Y are the same?

        int z0 = Math.round(nSlices/2);
        int nFramesToAverage = z0 - (int) Math.round(getMaxfromArray(offsetMFMarray)* axialSpacing/zStep);
        IJ.log("Number of frames to average: "+nFramesToAverage);
        IJ.log("z0: "+z0);
        int referenceFrameNumber = (int) Math.ceil((double) (nSplits*nSplits)/2);
        IJ.log("Reference frame: "+referenceFrameNumber);

        ImageStack imsAvg = new ImageStack(cropSizeX, cropSizeY, nSplits*nSplits);

        IJ.log("Cropping FOV...");
        int x,y,z;
        for (int j = 0; j < nSplits; j++){
            for (int i = 0; i < nSplits; i++){
                x = Math.round((width/nSplits)*i) + borderCrop;
                y = Math.round((height/nSplits)*j) + borderCrop;
                z = z0 + (int) Math.round(offsetMFMarray[i+j*nSplits]*axialSpacing/zStep) - nFramesToAverage;
//                IJ.log("x/y/z: "+x+"/"+y+"/"+z);
                ImageStack imsTemp = ims.crop(x, y, z, cropSizeX,cropSizeY, nFramesToAverage);
                imsAvg.setProcessor(calculateAverage(imsTemp), i+j*nSplits+1);
            }
        }

        ImagePlus impAvg = new ImagePlus("Avg stack", imsAvg);
        impAvg.setTitle(imp.getShortTitle()+" - Average");
        impAvg.show();

        ImageProcessor ipRef = imsAvg.getProcessor(referenceFrameNumber);
        GetShiftAndTiltRCCM RCCMcalculator = new GetShiftAndTiltRCCM(ipRef, nAngles, maxAngle);
        IJ.log("Calculating RCCM...");
        RCCMcalculator.computeRCCM(imsAvg);
        IJ.log("Extracting calibration parameters...");
        double[][] allShiftAndRot = RCCMcalculator.getShiftAndTiltfromRCCM();

        double[] shiftX = new double[nSplits*nSplits];
        double[] shiftY = new double[nSplits*nSplits];
        double[] theta = new double[nSplits*nSplits];

        for (int i = 0; i<nSplits*nSplits; i++){
            shiftX[i] = allShiftAndRot[0][i];
            shiftY[i] = allShiftAndRot[1][i];
            theta[i] = allShiftAndRot[2][i];
//            IJ.log("shift-X="+shiftX[i]+" shift-Y="+shiftY[i]+" tilt="+toDegrees(theta[i]));
        }

        ImagePlus impAvgCorr = applyCorrection(impAvg, shiftX, shiftY, theta);
        impAvgCorr.setTitle(imp.getShortTitle()+" - Average & Registered");
        impAvgCorr.show();

        IJ.log("Creating calibration table...");
        Map<String, double[]> data = new LinkedHashMap<>();
        data.put("X-shift (pixels)", shiftX);
        data.put("Y-shift (pixels)", shiftY);
        data.put("Theta (degrees)", theta); // in degrees !!!
        ResultsTable rt = dataMapToResultsTable(data);
        rt.show("Calibration-Table");
        try {
            saveNanoJTable(filePath+" - MFM_CalibrationTable.njt", rt);
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (applyCorr) {
            IJ.log("Applying correction to stack...");
            ImagePlus[] impCorrected = new ImagePlus[nSplits*nSplits];
            int maxOffset = getMaxfromArray(offsetMFMarray);
            double[] shiftXslice = new double[nSlices];
            double[] shiftYslice = new double[nSlices];
            double[] thetaSlice = new double[nSlices];
            for (int j = 0; j < nSplits; j++){
                for (int i = 0; i < nSplits; i++){

                    for (int k = 0; k < nSlices; k++){
                        shiftXslice[k] = shiftX[j*nSplits+i];
                        shiftYslice[k] = shiftY[j*nSplits+i];
                        thetaSlice[k] = theta[j*nSplits+i];
                    }
                    x = Math.round((width/nSplits)*i) + borderCrop;
                    y = Math.round((height/nSplits)*j) + borderCrop;
                    ImageStack imsTemp = ims.crop(x, y, 0, cropSizeX,cropSizeY, nSlices);
                    ImagePlus impTemp = new ImagePlus("Substack-"+i+"/"+j,imsTemp);
                    impCorrected[offsetMFMarray[j*nSplits+i]+maxOffset] = applyCorrection(impTemp, shiftXslice, shiftYslice, thetaSlice);
                }
            }
            Concatenator concatenator = new Concatenator();
            ImagePlus impStack = concatenator.concatenate(impCorrected, false);
            impStack.show();
            IJ.log("------------");
            IJ.log("All done.");

        }
    }

    // Helper function // TODO: could be migrated to core
    public static int getMaxfromArray(int[] array){
        int max = array[0];
        for (int i = 1; i < array.length; i++){
            if (max < array[i]) max = array[i];
        }
        return max;
    }
}

