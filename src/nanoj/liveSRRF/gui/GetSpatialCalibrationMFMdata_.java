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
import nanoj.liveSRRF.GetAxialPositionMFM;
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
    private final String MFMGetCalibVersion = "v0.3";

    @Override
    public void run(String s) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000 * imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep * 100d) / 100d;

        int[] offsetMFMarray = new int[] {2,3,4,-1,0,1,-4,-3,-2}; // TODO: add other cases where it's not 9 positions
        int nImageSplits = (int) Math.sqrt(offsetMFMarray.length);

        String[] splitLabels = new String[nImageSplits*nImageSplits];
        boolean[] splitDefaultChoice = new boolean[nImageSplits*nImageSplits];
        for (int i=0; i<nImageSplits*nImageSplits; i++){
            splitLabels[i] = Integer.toString(offsetMFMarray[i]);
            splitDefaultChoice[i] = true;
        }

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Spatial calibration from MFM data ("+MFMGetCalibVersion+")");
        gd.addMessage("-=-= MFM parameters =-=-\n");
        gd.addCheckboxGroup(nImageSplits, nImageSplits, splitLabels, splitDefaultChoice);
        gd.addNumericField("Axial spacing (in nm):", prefs.get("axialSpacing", 390), 1);
        gd.addNumericField("Border crop (in pixels):", prefs.get("borderCrop", 10), 0);
        gd.addMessage("-=-= Search parameters =-=-\n");
        gd.addNumericField("Number of angles:", prefs.get("nAngles", 20), 0);
        gd.addNumericField("Max. angle (in degrees):", prefs.get("maxAngle", 2), 1);
//        gd.addCheckbox("Apply correction", prefs.get("applyCorr", true));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int nROI = 0;
        boolean[] splitLocations = new boolean[nImageSplits*nImageSplits];
        for (int i=0; i<nImageSplits*nImageSplits; i++){
            splitLocations[i] = gd.getNextBoolean();
            if (splitLocations[i]) nROI++;
        }

        double[] nominalAxialPositions = new double[nROI];
        double[] chosenROIsLocations = new double[nROI];

        int id = 0;
        for (int i=0; i<nImageSplits*nImageSplits; i++){
            if (splitLocations[i]) {
                nominalAxialPositions[id] = offsetMFMarray[i];
                chosenROIsLocations[id] = i;
                id++;
            }
        }
        int[] sortedIndicesROI = getSortedIndices(nominalAxialPositions);

        float axialSpacing = (float) gd.getNextNumber();
        int borderCrop = (int) gd.getNextNumber();

        int nAngles = (int) gd.getNextNumber();
        float maxAngle = (float) gd.getNextNumber();
//        boolean applyCorr = gd.getNextBoolean();

//        prefs.set("nImageSplits", (float) nImageSplits);
        prefs.set("axialSpacing", axialSpacing);
        prefs.set("borderCrop", (float) borderCrop);
        prefs.set("nAngles", (float) nAngles);
        prefs.set("maxAngle", maxAngle);
//        prefs.set("applyCorr", applyCorr);

        String header = "Choose where to save Calibration table...";
        String filePath = getSavePath(header, imp.getTitle(), ".njt");
        if (filePath == null) return;
        filePath = filePath.replace(".njt", "");

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

        int cropSizeX = Math.round(width/nImageSplits) - 2*borderCrop; // TODO: assume that size in X and Y are the same?
        int cropSizeY = Math.round(height/nImageSplits) - 2*borderCrop; // TODO: assume that size in X and Y are the same?

        int z0 = Math.round(nSlices/2);
        int nFramesToAverage = z0 - (int) Math.round(getMaxfromArray(nominalAxialPositions) * axialSpacing/zStep); //one-sided: so 2x this
        IJ.log("Number of frames averaged: "+2*nFramesToAverage);
        IJ.log("z0: "+z0);
        int referenceFrameNumber = (int) Math.ceil((double) (nROI)/2)-1;
        IJ.log("Reference frame: "+referenceFrameNumber+" (position "+(int) chosenROIsLocations[referenceFrameNumber]+")");
        IJ.log("Number of angles: "+nAngles+" with Max. angle: "+maxAngle+" degrees ("+(2*maxAngle/(nAngles-1))+" degree(s) step)");

        ImageStack imsAvg = new ImageStack(cropSizeX, cropSizeY, nROI);

        IJ.log("------------");
        IJ.log("Cropping FOV...");
        int x,y,z;
        id = 0;
        for (int j = 0; j < nImageSplits; j++){
            for (int i = 0; i < nImageSplits; i++){
                if (splitLocations[j*nImageSplits+i]) {
                    x = Math.round((width / nImageSplits) * i) + borderCrop;
                    y = Math.round((height / nImageSplits) * j) + borderCrop;
                    z = z0 + (int) Math.round(offsetMFMarray[i + j * nImageSplits] * axialSpacing / zStep) - nFramesToAverage;
//                IJ.log("x/y/z: "+x+"/"+y+"/"+z);
                    ImageStack imsTemp = ims.crop(x, y, z, cropSizeX, cropSizeY, 2*nFramesToAverage);
                    imsAvg.setProcessor(calculateAverage(imsTemp), id+1);
                    id++;
                }
            }
        }

        ImagePlus impAvg = new ImagePlus("Avg stack", imsAvg);
        impAvg.setTitle(imp.getShortTitle()+" - Average");
        impAvg.show();

        ImageProcessor ipRef = imsAvg.getProcessor(referenceFrameNumber+1);
        GetShiftAndTiltRCCM RCCMcalculator = new GetShiftAndTiltRCCM(ipRef, nAngles, maxAngle);
        IJ.log("Calculating RCCM...");
        RCCMcalculator.computeRCCM(imsAvg);

//        ImageStack[] imsRCCMap = RCCMcalculator.imsRCCMap;
//        ImagePlus impRCCM1 = new ImagePlus("RCCM frame 1", imsRCCMap[0]);
//        impRCCM1.show();
//        ImagePlus impRCCM2 = new ImagePlus("RCCM frame 2", imsRCCMap[1]);
//        impRCCM2.show();

        IJ.log("Extracting calibration parameters...");
        double[][] allShiftAndRot = RCCMcalculator.getShiftAndTiltfromRCCM();

        double[] shiftX = new double[nROI];
        double[] shiftY = new double[nROI];
        double[] theta = new double[nROI];

        for (int i = 0; i<nROI; i++){
            shiftX[i] = allShiftAndRot[0][i];
            shiftY[i] = allShiftAndRot[1][i];
            theta[i] = allShiftAndRot[2][i];
//            IJ.log("shift-X="+shiftX[i]+" shift-Y="+shiftY[i]+" tilt="+toDegrees(theta[i]));
        }

        ImagePlus impAvgCorr = applyCorrection(impAvg, shiftX, shiftY, theta);
        impAvgCorr.setTitle(imp.getShortTitle()+" - Average & Registered");
        impAvgCorr.show();


        IJ.log("Applying correction to stack...");
        ImagePlus[] impCorrected = new ImagePlus[nROI];
        double[] shiftXslice = new double[nSlices];
        double[] shiftYslice = new double[nSlices];
        double[] thetaSlice = new double[nSlices];
        id = 0;
        for (int j = 0; j < nImageSplits; j++) {
            for (int i = 0; i < nImageSplits; i++) {
//                IJ.log("i/j: " + i + "/" + j + " (" + splitLocations[j * nImageSplits + i] + ")");
                if (splitLocations[j * nImageSplits + i]) {
                    for (int k = 0; k < nSlices; k++) {
                        shiftXslice[k] = shiftX[id];
                        shiftYslice[k] = shiftY[id];
                        thetaSlice[k] = theta[id];
                    }
                    x = Math.round((width / nImageSplits) * i) + borderCrop;
                    y = Math.round((height / nImageSplits) * j) + borderCrop;
                    ImageStack imsTemp = ims.crop(x, y, 0, cropSizeX, cropSizeY, nSlices);
                    ImagePlus impTemp = new ImagePlus("Substack-" + i + "/" + j, imsTemp);
                    impCorrected[sortedIndicesROI[id]] = applyCorrection(impTemp, shiftXslice, shiftYslice, thetaSlice);
                    id++;
                }
            }
        }

//        IJ.log("impCorr size: "+impCorrected.length+" id="+id);

        IJ.log("Getting axial positions...");
        ImageStack imsRef = impCorrected[referenceFrameNumber + 1].getImageStack();
        imsRef = imsRef.crop(0, 0, z0 - nFramesToAverage/2, cropSizeX, cropSizeY, nFramesToAverage); // TODO: which chunk is best? currently divided by 2
        GetAxialPositionMFM getAxialPositionMFM = new GetAxialPositionMFM(imsRef, (float) zStep, nSlices);
        double[] axialPositions = new double[nROI];
        for (int i = 0; i < nROI; i++) {
            axialPositions[i] = getAxialPositionMFM.computeZcorrelation(impCorrected[i].getImageStack());
        }

        double meanOffset = 0;
        double axialPositionRef = axialPositions[referenceFrameNumber];
        for (int i = 0; i < nROI; i++) {
            axialPositions[i] -= axialPositionRef;
            if (i!=referenceFrameNumber) meanOffset += axialPositions[i]/nominalAxialPositions[sortedIndicesROI[i]];
        }
        meanOffset /= nROI-1;

        IJ.log("Creating calibration table...");
        Map<String, double[]> data = new LinkedHashMap<>();
        data.put("X-shift (pixels)", sortArray(shiftX, sortedIndicesROI));
        data.put("Y-shift (pixels)", sortArray(shiftY, sortedIndicesROI));
        data.put("Theta (degrees)", sortArray(theta, sortedIndicesROI)); // in degrees !!!
        data.put("Axial positions", axialPositions);
        data.put("Nominal positions", sortArray(nominalAxialPositions, sortedIndicesROI));
        data.put("ROI #", sortArray(chosenROIsLocations, sortedIndicesROI));
        ResultsTable rt = dataMapToResultsTable(data);
        rt.show("Calibration-Table");
        try {
            saveNanoJTable(filePath+" - MFM_CalibrationTable.njt", rt);
        } catch (IOException e) {
            e.printStackTrace();
        }

        Concatenator concatenator = new Concatenator();
        ImagePlus impStack = concatenator.concatenate(impCorrected, false);
        impStack.setTitle(imp.getShortTitle()+" - Calibrated");
        impStack.show();

        IJ.log("------------");
        IJ.log("Mean offset: "+meanOffset+" nm");
        IJ.log("------------");
        IJ.log("All done.");

    }

    // Helper function // TODO: could be migrated to core
    public static double getMaxfromArray(double[] array){
        double max = array[0];
        for (int i = 1; i < array.length; i++){
            if (max < array[i]) max = array[i];
        }
        return max;
    }

    // Helper function // TODO: could be migrated to core
    public static int[] getSortedIndices(double[] array){
        int[] indices = new int[array.length];
        for (int i=0; i<array.length; i++){
            for (double v : array) {
                if (v < array[i]) indices[i]++;
            }
        }
        return indices;
    }

    public static double[] sortArray(double[] array, int[] indices){
        double[] sortedArray = new double[array.length];

        for (int i=0; i<array.length; i++){
            sortedArray[i] = array[indices[i]];
        }
        return sortedArray;
    }

}

