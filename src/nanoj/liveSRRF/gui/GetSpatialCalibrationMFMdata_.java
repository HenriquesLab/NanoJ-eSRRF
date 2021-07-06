package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.io.FileInfo;
import ij.measure.ResultsTable;
import ij.plugin.HyperStackConverter;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.array.ArrayMath;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.GetAxialPositionMFM;
import nanoj.liveSRRF.MFMCalibration;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import static nanoj.core.java.imagej.FilesAndFoldersTools.getSavePath;
import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.core.java.io.SaveNanoJTable.saveNanoJTable;
import static nanoj.liveSRRF.LinearRegressions.linearRegressionLeastSquareNoOffset;
import static nanoj.liveSRRF.StackProjections.calculateAverage;
import static nanoj.liveSRRF.StackProjections.calculateMIP;

public class GetSpatialCalibrationMFMdata_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    private final String MFMGetCalibVersion = "v0.7";

    @Override
    public void run(String s) {

        // ---- Getting input data ----
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000 * imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep * 100d) / 100d;

        // ---- Building the GUI ----
        int[] offsetMFMarray = new int[] {2,3,4,-1,0,1,-4,-3,-2}; // TODO: add other cases where it's not 9 positions? or get the positions automatically?
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
        gd.addMessage("-=-= Intensity rescaling =-=-\n");
        gd.addNumericField("Background level:", prefs.get("bgLevel", 100), 2);
        gd.addNumericField("Percentile", prefs.get("pcTile", 1),2);
        gd.addMessage("-=-= Display options =-=-\n");
        gd.addCheckbox("Show RCCMs", prefs.get("showRCCMs", false));
        gd.addCheckbox("Show intermediate results", prefs.get("showAllResults", false));
        gd.addCheckbox("Show plots", prefs.get("showPlots", false));


        gd.showDialog();
        if (gd.wasCanceled()) return;

        // ---- Getting info from GUI ----
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
        float bgLevel = (float) gd.getNextNumber();
        float pcTile = (float) gd.getNextNumber();
        boolean showRCCMs = gd.getNextBoolean();
        boolean showAllResults = gd.getNextBoolean();
        boolean showPlots = gd.getNextBoolean();

        // ---- Saving GUI input for next use ----
        prefs.set("axialSpacing", axialSpacing);
        prefs.set("borderCrop", (float) borderCrop);
        prefs.set("nAngles", (float) nAngles);
        prefs.set("maxAngle", maxAngle);
        prefs.set("bgLevel", bgLevel);
        prefs.set("pcTile", pcTile);
        prefs.set("showRCCMs", showRCCMs);
        prefs.set("showAllResults", showAllResults);
        prefs.set("showPlots", showPlots);

        // ---- Get path for NanoJ table save ----
        String header = "Choose where to save Calibration table...";
        String filePath = getSavePath(header, imp.getTitle(), ".njt");
        if (filePath == null) return;
        filePath = filePath.replace(".njt", "");

        // ---- Initialize the variables and the Log ----
        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("MFM Calibration "+MFMGetCalibVersion);
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("Axial scan step size (in nm): "+zStep);

        ImageStack ims = imp.duplicate().getImageStack().convertToFloat();
        int nSlices = ims.getSize();
        int width = ims.getWidth();
        int height = ims.getHeight();

        // ---- Removing background level from image stack ----
        for (int i=0; i<ims.getSize(); i++){
            ims.getProcessor(i+1).add( -bgLevel);
        }

        // ---- Initialize the variables and the Log ----
        int cropSizeX = width/nImageSplits - 2*borderCrop; // int division rounds down automatically, ignoring the remainder
        int cropSizeY = height/nImageSplits - 2*borderCrop;
//IJ.log("CropSizeX: "+cropSizeX);

        int z0 = Math.round((float) nSlices/2);
        int nFramesToAverage = z0 - (int) Math.ceil(getMaxfromArray(nominalAxialPositions) * axialSpacing / zStep + 1); //one-sided: so 2x this
        IJ.log("Number of frames averaged: "+2*nFramesToAverage);
        IJ.log("z0: "+z0);
        int referenceFrameNumber = (int) Math.ceil((double) (nROI)/2)-1; // this is zero-based
        IJ.log("Reference frame: "+referenceFrameNumber+" (position "+(int) chosenROIsLocations[referenceFrameNumber]+")");
        IJ.log("Number of angles: "+nAngles+" with Max. angle: "+maxAngle+" degrees ("+(2*maxAngle/(nAngles-1))+" degree(s) step)");

        ImageStack imsAvg = new ImageStack(cropSizeX, cropSizeY, nROI);
        ImageStack imsMIP = new ImageStack(cropSizeX, cropSizeY, nROI);

        // ---- Cropping loop based on nominal xyz positions ----
        IJ.log("------------");
        IJ.log("Cropping FOV...");
        int x,y,z;
        id = 0;
        for (int j = 0; j < nImageSplits; j++){
            for (int i = 0; i < nImageSplits; i++){
                if (splitLocations[j*nImageSplits+i]) {
                    x = Math.round((float) width / nImageSplits * i) + borderCrop;
                    y = Math.round((float) height / nImageSplits * j) + borderCrop;
                    z = z0 + Math.round(offsetMFMarray[i + j * nImageSplits] * axialSpacing / (float) zStep) - nFramesToAverage;
//                IJ.log("x/y/z: "+x+"/"+y+"/"+z);
                    ImageStack imsTemp = ims.crop(x, y, z, cropSizeX, cropSizeY, 2*nFramesToAverage);
                    imsAvg.setProcessor(calculateAverage(imsTemp), id+1);
                    imsMIP.setProcessor(calculateMIP(imsTemp), id+1);
                    id++;
                }
            }
        }

        // ---- Showing average data stacks prior to registration ----
        if (showAllResults){
            ImagePlus impAvg = new ImagePlus(imp.getShortTitle()+" - Average", imsAvg);
            impAvg.show();
        }

        // ---- Calculation xy shift and theta rotation ----
        ImageProcessor ipRef = imsMIP.getProcessor(referenceFrameNumber+1);
        MFMCalibration RCCMcalculator = new MFMCalibration();

        IJ.log("Calculating RCCM...");
        RCCMcalculator.computeRCCM(ipRef, imsMIP, nAngles, maxAngle); // angles are fed in degrees


        // -----------
        if (showRCCMs) {
            ImageStack[] imsRCCM = RCCMcalculator.imsRCCMap;
            for (int i = 0; i < imsRCCM.length; i++) {
                ImagePlus impRCCMTemp = new ImagePlus("RCCM" + i, imsRCCM[i]);
                impRCCMTemp.show();
            }
        }
        // -----------

        IJ.log("Extracting calibration parameters...");
        double[][] allShiftAndRot = RCCMcalculator.getShiftAndTiltfromRCCM(showPlots);
        double[] shiftX = new double[nROI];
        double[] shiftY = new double[nROI];
        double[] theta = new double[nROI];

        for (int i = 0; i<nROI; i++){
            shiftX[i] = allShiftAndRot[0][i];
            shiftY[i] = allShiftAndRot[1][i];
            theta[i] = allShiftAndRot[2][i];
//            IJ.log("shift-X="+shiftX[i]+" shift-Y="+shiftY[i]+" tilt="+toDegrees(theta[i]));
        }

        // ---- Applying corrections to average data stack ----
        ImageStack[] imsCorr = RCCMcalculator.applyMFMCorrection(imsAvg, shiftX, shiftY, theta, null, null);

        if (showAllResults) {
            ImagePlus impAvgCorr = new ImagePlus(imp.getShortTitle() + " - Average & Registered", imsCorr[0]);
            impAvgCorr.show();
        }

        if (showAllResults) {
            ImagePlus impAvgCorrCropped = new ImagePlus(imp.getShortTitle() + " - Average & Registered & Cropped", imsCorr[1]);
            impAvgCorrCropped.show();
        }

        // This uses the Cropped and Registered imsAvg
        double[] intensityScalingCoeffs = getLinearRegressionParameters(imsCorr[1], referenceFrameNumber, (100-pcTile)/100);

        if (showAllResults) {
            ImageStack imsAvgCorrCropRescaled = imsCorr[1].duplicate();
            for (int i = 0; i < nROI; i++) {
                imsAvgCorrCropRescaled.getProcessor(i+1).multiply(1/intensityScalingCoeffs[i]);
            }
            ImagePlus impAvgCorrCropRescaled = new ImagePlus(imp.getShortTitle() + " - Average & Registered & Cropped & Rescaled", imsAvgCorrCropRescaled);
            impAvgCorrCropRescaled.show();
        }

        // ---- Applying correction to the whole stack stack ----
        IJ.log("Applying corrections to stack...");
        ImageStack[] imsCorrected = new ImageStack[nROI];
        ImageStack imsCorrectedUberStack = new ImageStack();

        double[] shiftXslice = new double[nSlices];
        double[] shiftYslice = new double[nSlices];
        double[] thetaSlice = new double[nSlices];
        double[] coeffSlice = new double[nSlices];

        id = 0;
        for (int j = 0; j < nImageSplits; j++) {
            for (int i = 0; i < nImageSplits; i++) {
//                IJ.log("i/j: " + i + "/" + j + " (" + splitLocations[j * nImageSplits + i] + ")");
                if (splitLocations[j * nImageSplits + i]) {
                    for (int k = 0; k < nSlices; k++) {
                        shiftXslice[k] = shiftX[id];
                        shiftYslice[k] = shiftY[id];
                        thetaSlice[k] = theta[id];
                        coeffSlice[k] = intensityScalingCoeffs[id];
                    }
                    x = Math.round((float) width / nImageSplits * i) + borderCrop;
                    y = Math.round((float) height / nImageSplits * j) + borderCrop;
                    ImageStack imsTemp = ims.crop(x, y, 0, cropSizeX, cropSizeY, nSlices);
                    imsCorrected[sortedIndicesROI[id]] = RCCMcalculator.applyMFMCorrection(imsTemp, shiftXslice, shiftYslice, thetaSlice, coeffSlice, null)[0];
                    id++;
                }
            }
        }

        // ---- Obtaining the relative axial positions ----
        IJ.log("Getting axial positions...");
        ImageStack imsRef = imsCorrected[referenceFrameNumber]; // referenceFrameNumber is zero-based so all good here
        imsRef = imsRef.crop(0, 0, z0 - nFramesToAverage/2, cropSizeX, cropSizeY, nFramesToAverage); // TODO: which chunk is best? currently divided by 2
        GetAxialPositionMFM getAxialPositionMFM = new GetAxialPositionMFM(imsRef, (float) zStep, nSlices);
        double[] axialPositions = new double[nROI];
        for (int i = 0; i < nROI; i++) {
            axialPositions[i] = getAxialPositionMFM.computeZcorrelation(imsCorrected[i], showPlots);
        }

        // ---- Calculating the mean offset ----
        double axialPositionRef = axialPositions[referenceFrameNumber];
        for (int i = 0; i < nROI; i++) {
            axialPositions[i] -= axialPositionRef;
        }
        double meanOffset = linearRegressionLeastSquareNoOffset(sortArray(nominalAxialPositions, sortedIndicesROI), axialPositions);


        // ---- Creating the NanoJ table and save it ----
        int digitsToRound = 2; // should be fine
        IJ.log("Creating calibration table...");
        Map<String, double[]> data = new LinkedHashMap<>();
        data.put("X-shift (pixels)", sortArray(roundToNdigits(shiftX, digitsToRound), sortedIndicesROI));
        data.put("Y-shift (pixels)", sortArray(roundToNdigits(shiftY, digitsToRound), sortedIndicesROI));
        data.put("Theta (degrees)", sortArray(roundToNdigits(theta, digitsToRound), sortedIndicesROI)); // in degrees !!!
        data.put("Axial positions", roundToNdigits(axialPositions, digitsToRound));
        data.put("Nominal positions", sortArray(roundToNdigits(nominalAxialPositions, digitsToRound), sortedIndicesROI));
        data.put("ROI #", sortArray(roundToNdigits(chosenROIsLocations, digitsToRound), sortedIndicesROI));
        data.put("Intensity scaling", sortArray(roundToNdigits(intensityScalingCoeffs, digitsToRound), sortedIndicesROI));
        double[] bgArray = new double[nROI];
        for (int i = 0; i < nROI; i++) {
            bgArray[i] = bgLevel;
        }
        data.put("Background level (ADC)", bgArray);

        ResultsTable rt = dataMapToResultsTable(data);
        rt.show("Calibration-Table");
        try {
            saveNanoJTable(filePath+" - MFM_CalibrationTable.njt", rt);
        } catch (IOException e) {
            e.printStackTrace();
        }

        // ---- Reshape the stack into Hyperstack ----
        for (int i = 0; i < nROI; i++) {
            for (int k = 0; k < nSlices; k++) {
                imsCorrectedUberStack.addSlice(imsCorrected[i].getProcessor(k+1));
            }
        }

        ImagePlus impStack = new ImagePlus(imp.getShortTitle()+" - Calibrated", imsCorrectedUberStack);
        impStack.setTitle(imp.getShortTitle()+" - Calibrated");
        impStack = HyperStackConverter.toHyperStack(impStack, 1, nROI, nSlices, "xyctz", "Color");
        impStack.show();


        IJ.log("------------");
        IJ.log("Offset: "+roundToNdigits(new double[]{meanOffset}, digitsToRound)[0]+" nm");
        IJ.log("------------");
        IJ.log("All done.");

    }


    // ----------------------------------------- FUNCTIONS -----------------------------------------
    // Helper function // TODO: could be migrated to core
    public static double getMaxfromArray(double[] array){
        double max = array[0];
        for (int i = 1; i < array.length; i++){
            if (max < array[i]) max = array[i];
        }
        return max;
    }

    // Helper function // TODO: could be migrated to core, assumes that there is an order to be found and no repeats
    public static int[] getSortedIndices(double[] array){
        int[] indices = new int[array.length];
        for (int i=0; i<array.length; i++){
            for (double v : array) {
                if (v < array[i]) indices[i]++;
            }
        }
        return indices;
    }

    // Helper function // TODO: could be migrated to core
    public static double[] sortArray(double[] array, int[] indices){
        double[] sortedArray = new double[array.length];

        for (int i=0; i<array.length; i++){
            sortedArray[i] = array[indices[i]];
        }
        return sortedArray;
    }

    // Helper function // TODO: could be migrated to core
    public static double[] getLinearRegressionParameters(ImageStack ims, int refSlice, float pcTile){

        ImageStack[] imsMasks = getMasksFromIntensityPercentile(ims.duplicate(), pcTile);
        FloatProcessor fp = imsMasks[1].getProcessor(1).convertToFloatProcessor();
        float[] pixelMask = (float[]) fp.getPixels(); // Using getVOxels does not work directly on the ims if no imp is connected to it
//        IJ.log("Sum of pixels: "+ArrayMath.sum(pixelMask));

        // This uses a least square-based linear regression assuming intercept of 0
        float[] arrayRef = new float[ims.getWidth()*ims.getHeight()];
        float[] arrayTemp = new float[ims.getWidth()*ims.getHeight()];
        arrayRef = ims.getVoxels(0,0, refSlice, ims.getWidth(), ims.getHeight(), 1, arrayRef);

        float XX = 0;
        for (int i = 0; i < ims.getWidth()*ims.getHeight(); i++) {
            XX += arrayRef[i]*arrayRef[i]*pixelMask[i];


        }
//        IJ.log("XX="+XX);

        double[] coeffs = new double[ims.getSize()];
        for (int s = 0; s < ims.getSize(); s++) {
            arrayTemp = ims.getVoxels(0,0, s, ims.getWidth(), ims.getHeight(), 1, arrayTemp);

            for (int i = 0; i < ims.getWidth()*ims.getHeight(); i++) {
                coeffs[s] += arrayRef[i]*arrayTemp[i]*pixelMask[i];
            }
            coeffs[s] /= XX;
//            IJ.log("Coeffs: "+coeffs[i]);
        }
        return coeffs;
    }

    // Helper function ---
    public static ImageStack[] getMasksFromIntensityPercentile(ImageStack ims, float percentTile){
// This code returns an array of ImageStack, the first element is an image stack which contains all the individual Masks
// Each mask is calculated based on their individual percentile
// The second element returned is a stack with the mask representing the overlap of all the masks
        ImageStack imsMasks = new ImageStack(ims.getWidth(), ims.getHeight(), ims.getSize());
        float[] sortedArrayTemp;
        float[] arrayOverlappingMask = new float[ims.getWidth()*ims.getHeight()];
        float valueAtpcTile;
        float vTemp;

        if (percentTile < 0) percentTile = 0;
        if (percentTile > 1) percentTile = 1;

        for (int s = 0; s < ims.getSize(); s++) {
            float[] arrayTemp = new float[ims.getWidth()*ims.getHeight()];
            arrayTemp = ims.duplicate().getVoxels(0,0, s, ims.getWidth(), ims.getHeight(), 1, arrayTemp);
            sortedArrayTemp = arrayTemp.clone();
            Arrays.sort(sortedArrayTemp);
            valueAtpcTile = sortedArrayTemp[(int) ((arrayTemp.length-1) * percentTile)];
//            IJ.log("Value at "+(100*percentTile)+"th percentile: "+valueAtpcTile);
            for (int y = 0; y < ims.getHeight(); y++) {
                for (int x = 0; x < ims.getWidth(); x++) {
                    if (arrayTemp[x + y*ims.getWidth()] > valueAtpcTile) vTemp = 1;
                    else vTemp = 0;
//                    imsMasks.setVoxel(x, y, s, vTemp);
                    arrayTemp[x + y*ims.getWidth()] = vTemp;
                    if (s == 0) arrayOverlappingMask[x + y*ims.getWidth()] = vTemp;
                    else arrayOverlappingMask[x + y*ims.getWidth()] *= vTemp;
                }
            }
            imsMasks.setProcessor(new FloatProcessor(ims.getWidth(), ims.getHeight(), arrayTemp), s+1);
        }

        IJ.log("Number of pixels in mask: "+ ArrayMath.sum(arrayOverlappingMask));
        ImageStack imsOverlappingMask = new ImageStack(ims.getWidth(), ims.getHeight(), 1);
        imsOverlappingMask.setProcessor(new FloatProcessor(ims.getWidth(), ims.getHeight(), arrayOverlappingMask), 1);
        return new ImageStack[]{imsMasks, imsOverlappingMask};
    }

    // Helper function --- this rounds an array to "digits" digits
    public static double[] roundToNdigits(double[] array, int digits){

        double R = Math.pow(10, digits);
        double[] roundedArray = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            roundedArray[i] = Math.round(array[i]*R)/R;
        }
        return roundedArray;
    }


}

