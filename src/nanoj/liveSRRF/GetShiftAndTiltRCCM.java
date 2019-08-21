package nanoj.liveSRRF;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui.ApplyDriftCorrection;
import nanoj.core.java.image.transform.TranslateOrRotateImage;
import nanoj.core.java.tools.Log;

import java.awt.*;

import static java.lang.Math.*;
import static nanoj.core.java.array.ArrayInitialization.initializeDoubleAndGrowthFill;
import static nanoj.core.java.array.ArrayMath.getMaxValue;
import static nanoj.core.java.array.ArrayMath.getMinValue;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.array.ArrayInitialization.initializeFloatAndGrowthFill;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getMaxFindByOptimization;
import static nanoj.liveSRRF.GetAxialPositionMFM.normalizeArray;

public class GetShiftAndTiltRCCM {
    private static final TranslateOrRotateImage translateOrRotateImage = new TranslateOrRotateImage();

    private ImageProcessor ipRef;
    private float[] angleArray;
    private double[] angleArrayDouble;
    private float angleStep;
    private static final Log log = new Log();
    //    private int radiusX, radiusY; // TODO: include radius X and Y to speed things up?
    public ImageStack[] imsRCCMap;
    Plot thetaFitPlot;

    // Initialisation
    public GetShiftAndTiltRCCM(ImageProcessor ipRef, int nAngles, float maxAngle){
        this.ipRef = ipRef.duplicate();
        this.angleStep = 2*(float) toRadians(maxAngle)/(nAngles-1);
        this.angleArray = initializeFloatAndGrowthFill(nAngles, -(float) toRadians(maxAngle), angleStep);
        this.angleArrayDouble = initializeDoubleAndGrowthFill(nAngles, -toRadians(maxAngle), angleStep);
        this.thetaFitPlot = new Plot("Theta fit plot", "Angle (radians)", "Intensity (AU)");
    }

    // This method calculates the RCCM ImageStack array ----------------------------------------------------------------
    public void computeRCCM(ImageStack ims) {
        int nSlices = ims.getSize();
        imsRCCMap = new ImageStack[nSlices];

        // Rotate and then calculate Cross Correlation
        for (int s = 1; s <= nSlices; s++) {

            log.progress(s, nSlices);

            // build rotated stack
            ImageStack imsRotated = new ImageStack(ipRef.getWidth(), ipRef.getHeight());
            for (int n = 0; n < angleArray.length; n++) {
                imsRotated.addSlice(ims.getProcessor(s));
            }

            // rotate
            imsRotated = translateOrRotateImage.rotate(imsRotated, angleArray);
            imsRCCMap[s - 1] = calculateCrossCorrelationMap(ipRef, imsRotated, false);
        }
    }

    // This method calculates the shift and tilt from the RCCM ImageStack array ----------------------------------------
    public double[][] getShiftAndTiltfromRCCM(){

        int nROI = imsRCCMap.length;
        Fit1DGaussian fitting;
        float cropLevel = 0.2f;

        double[] shiftX = new double[nROI];
        double[] shiftY = new double[nROI];
        double[] theta = new double[nROI];

        int nSlices;

        for(int s=0; s<nROI; s++) {
            ImageStack imsRCCMapSlice = imsRCCMap[s];

            // calculate position and value of maximum pixel TODO: avoid spurious pixels that would skew this
            nSlices = imsRCCMapSlice.getSize();
            float[][] shiftXY = new float[2][nSlices];
            float[] maxValue = new float[nSlices];

            for (int a = 0; a<nSlices; a++) {
                FloatProcessor fpRCCMslice = imsRCCMapSlice.getProcessor(a+1).convertToFloatProcessor();
                float[] resultsOfOpimisation = getMaxFindByOptimization(fpRCCMslice);
                shiftXY[0][a] = resultsOfOpimisation[0];
                shiftXY[1][a] = resultsOfOpimisation[1];
                maxValue[a] = resultsOfOpimisation[2];
            }

            fitting = new Fit1DGaussian(normalizeArray(maxValue));
            fitting.cropDataArray(cropLevel);
            float[] fitResults = fitting.calculate();
            double[][] modelArray = fitting.fittedCurve();
//            IJ.log("--- Fit results ---");
//            IJ.log("Amp: "+fitResults[0]);
//            IJ.log("x0: "+(fitResults[1]));
//            IJ.log("Sigma: "+(fitResults[2]));
//            IJ.log("BG: "+fitResults[3]);

            theta[s] = toDegrees(linearInterpolation(angleArray, fitResults[1]));
//            shiftX[s] = (imsRCCMapSlice.getWidth()-1)/2 - shiftXY[0][Math.round(fitResults[1])];
//            shiftY[s] = (imsRCCMapSlice.getHeight()-1)/2 - shiftXY[1][Math.round(fitResults[1])];
            shiftX[s] = (imsRCCMapSlice.getWidth()-1)/2 - linearInterpolation(shiftXY[0], fitResults[1]); // TODO: this assumes that shiftXY is relatively continuous wrt to theta
            shiftY[s] = (imsRCCMapSlice.getHeight()-1)/2 - linearInterpolation(shiftXY[1], fitResults[1]);

            thetaFitPlot.setColor(Color.black);
            thetaFitPlot.add("line", angleArrayDouble, modelArray[0]);
            thetaFitPlot.setColor(Color.red);
            thetaFitPlot.add("line", angleArrayDouble, modelArray[1]);
//        defocusPlot.add("line", zPosArray, normalizeArray(zCorrArray));
            thetaFitPlot.show();
        }

        return new double[][] {shiftX, shiftY, theta};
    }

    public static ImagePlus[] applyCorrection(ImagePlus imp, double[] shiftX, double[] shiftY, double[] theta, double[] intCoeffs){

        float[] floatTheta = new float[theta.length];
        for (int i = 0 ; i < theta.length; i++) {
            floatTheta[i] = (float) toRadians(theta[i]);
        }

        ApplyDriftCorrection adc = new ApplyDriftCorrection();
        ImageStack ims = imp.getImageStack();

        if (intCoeffs != null){
            for (int i = 0; i < ims.getSize(); i++) {
                ims.getProcessor(i+1).multiply( 1/intCoeffs[i]);
            }
        }

        ImageStack imsRotated = translateOrRotateImage.rotate(ims, floatTheta);

        ImagePlus impRot = new ImagePlus(imp.getShortTitle(), imsRotated);
        ImagePlus impAvgCorrected = adc.applyDriftCorrection(impRot, shiftX, shiftY);

        int x0 = (int) Math.max(Math.ceil(getMaxValue(shiftX)[1]), -Math.ceil(getMinValue(shiftX)[1]));
        int y0 = (int) Math.max(Math.ceil(getMaxValue(shiftY)[1]), -Math.ceil(getMinValue(shiftY)[1]));
        ImageStack imsAvgCorrectedCropped = impAvgCorrected.getStack().crop(x0, y0, 0, imp.getWidth()-2*x0, imp.getHeight()-2*y0, imp.getStackSize());
        ImagePlus[] impCombo = new ImagePlus[2];
        impCombo[0] = impAvgCorrected;
        impCombo[1] = new ImagePlus("Corrected crop", imsAvgCorrectedCropped);

        return impCombo;
    }

    public static float linearInterpolation(float[] array, float x){

        if (x<0) return array[0];
        else if (x>array.length-1) return array[array.length-1];
        else {
            int x0 = (int) Math.floor((double) x);
            return ((1-x+x0)*array[x0]+(x-x0)*array[x0+1]);}

    }


}