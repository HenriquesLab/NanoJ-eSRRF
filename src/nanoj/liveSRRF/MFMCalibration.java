package nanoj.liveSRRF;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui.ApplyDriftCorrection;
import nanoj.core.java.image.transform.TranslateOrRotateImage;
import nanoj.core.java.threading.NanoJThreadExecutor;
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

public class MFMCalibration {

    private double[] angleArray;
    private static final Log log = new Log();
    //    private int radiusX, radiusY; // TODO: include radius X and Y to speed things up?
    public ImageStack[] imsRCCMap;

    // Initialisation
    public MFMCalibration(){
        // Keep calm and carry on
    }

    // This method calculates the RCCM ImageStack array ----------------------------------------------------------------
    public void computeRCCM(ImageProcessor ipRef, ImageStack ims, int nAngles, float maxAngle) {

        int nROI = ims.getSize();
        imsRCCMap = new ImageStack[nROI];
        if (nAngles <= 0 | maxAngle == 0){
            nAngles = 1;
            angleArray = new double[]{0};
        }
        else {
            double angleStep = 2 * maxAngle / (nAngles - 1); // in degrees
            angleArray = initializeDoubleAndGrowthFill(nAngles, -maxAngle, angleStep); // in degrees
        }

        // Rotate and then calculate Cross Correlation
        for (int s = 1; s <= nROI; s++) {

            log.progress(s, nROI);

            // build rotated stack
            ImageStack imsRotated = new ImageStack(ipRef.getWidth(), ipRef.getHeight());
            for (int n = 0; n < nAngles; n++) {
                imsRotated.addSlice(ims.getProcessor(s));
            }

            // rotate
            NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
            NTE.showProgress = false;
            for (int n = 1; n <= nAngles; n++) {
                ThreadedRotate r = new ThreadedRotate(imsRotated, n, angleArray[n - 1]);
                NTE.execute(r);
            }
            NTE.finish();

            // --- Using aparapi ---
//            imsRotated = translateOrRotateImage.rotate(imsRotated, angleArray);

            // --- Manual rotation ---
//            for (int n = 0; n < nAngles; n++) {
//                ImageProcessor ip = imsRotated.getProcessor(n+1).duplicate();
//                ip.setInterpolationMethod(ip.BICUBIC);
//                ip.rotate(toDegrees(angleArrayDouble[n]));
//                imsRotated.setProcessor(ip, n+1);
//            }

            imsRCCMap[s - 1] = calculateCrossCorrelationMap(ipRef, imsRotated,false);
        }
    }

    // This method calculates the shift and tilt from the RCCM ImageStack array ----------------------------------------
    public double[][] getShiftAndTiltfromRCCM(boolean showPlot){

        int nROI = imsRCCMap.length;
        Fit1DGaussian fitting;
        float cropLevel = 0.7f; // TODO: make this user-set?

        double[] shiftX = new double[nROI];
        double[] shiftY = new double[nROI];
        double[] theta = new double[nROI];

        int nAngles;
        Plot thetaFitPlot = new Plot("Theta fit plot", "Angle (degrees)", "Intensity (AU)");

        for(int s=0; s<nROI; s++) {
            ImageStack imsRCCMapSlice = imsRCCMap[s];

            // calculate position and value of maximum pixel
            nAngles = imsRCCMapSlice.getSize();
            float[][] shiftXY = new float[2][nAngles];
            float[] maxValue = new float[nAngles];

            for (int a = 0; a<nAngles; a++) {
                FloatProcessor fpRCCMslice = imsRCCMapSlice.getProcessor(a+1).convertToFloatProcessor();
                float[] resultsOfOpimisation = getMaxFindByOptimization(fpRCCMslice);
                shiftXY[0][a] = resultsOfOpimisation[0];
                shiftXY[1][a] = resultsOfOpimisation[1];
                maxValue[a] = resultsOfOpimisation[2];
            }

            if (nAngles > 1) {
                fitting = new Fit1DGaussian(normalizeArray(maxValue));
                fitting.cropDataArray(cropLevel);
                float[] fitResults = fitting.calculate();
                double[][] modelArray = fitting.fittedCurve();
//            IJ.log("--- Fit results ---");
//            IJ.log("Amp: "+fitResults[0]);
//            IJ.log("x0: "+(fitResults[1]));
//            IJ.log("Sigma: "+(fitResults[2]));
//            IJ.log("BG: "+fitResults[3]);

                theta[s] = linearInterpolation(angleArray, fitResults[1]); // convert it back to degrees
//            shiftX[s] = (imsRCCMapSlice.getWidth()-1)/2 - shiftXY[0][Math.round(fitResults[1])];
//            shiftY[s] = (imsRCCMapSlice.getHeight()-1)/2 - shiftXY[1][Math.round(fitResults[1])];
                shiftX[s] = (imsRCCMapSlice.getWidth() - 1) / 2 - linearInterpolation(shiftXY[0], fitResults[1]); // This assumes that shiftXY is relatively continuous wrt to theta
                shiftY[s] = (imsRCCMapSlice.getHeight() - 1) / 2 - linearInterpolation(shiftXY[1], fitResults[1]);

                if (showPlot) {
                    thetaFitPlot.setColor(Color.black);
                    thetaFitPlot.add("line", angleArray, modelArray[0]);
                    thetaFitPlot.setColor(Color.red);
                    thetaFitPlot.add("line", angleArray, modelArray[1]);
//        defocusPlot.add("line", zPosArray, normalizeArray(zCorrArray));
                    thetaFitPlot.show();
                }
            }
            else {
                theta[s] = angleArray[0]; // don't fit if there's only one angle in the array
                shiftX[s] = (imsRCCMapSlice.getWidth() - 1) / 2 - shiftXY[0][0];
                shiftY[s] = (imsRCCMapSlice.getHeight() - 1) / 2 - shiftXY[1][0];
            }
        }

        return new double[][] {shiftX, shiftY, theta};
    }

    // Threaded version which does not use Aparapi
    public ImageStack[] applyMFMCorrection(ImageStack ims, double[] shiftX, double[] shiftY, double[] theta, double[] intCoeffs, double[] bgLevels){

        ImageStack imsCorr = ims.duplicate();
        int nSlices = ims.getSize();

        if (intCoeffs != null){
//            IJ.log("C");
            for (int i = 0; i < nSlices; i++) {
                if (bgLevels != null) {
//                    IJ.log("BG");
                    imsCorr.getProcessor(i + 1).add(-bgLevels[i]);
                }
                imsCorr.getProcessor(i+1).multiply( 1.0d/intCoeffs[i]); // rescaling should only be done on BG-removed dataset!!
            }
        }

        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
        NTE.showProgress = false;
        for (int n = 1; n <= nSlices; n++) {
            log.status("Translating frame " + n + "/" + nSlices);
            log.progress(n, nSlices);
            ThreadedRotateTranslate rt = new ThreadedRotateTranslate(imsCorr, n, shiftX[n - 1], shiftY[n - 1], theta[n - 1]);
            NTE.execute(rt);
        }
        NTE.finish();

        int x0 = (int) Math.max(Math.ceil(getMaxValue(shiftX)[1]), -Math.ceil(getMinValue(shiftX)[1]));
        int y0 = (int) Math.max(Math.ceil(getMaxValue(shiftY)[1]), -Math.ceil(getMinValue(shiftY)[1]));
        ImageStack imsCorrCropped = imsCorr.crop(x0, y0, 0, ims.getWidth()-2*x0, ims.getHeight()-2*y0, nSlices);

        return new ImageStack[]{imsCorr, imsCorrCropped};
    }


    // -------------- Santa's little helper --------------
    // TODO: could be moved to NanoJ-Core
    public static float linearInterpolation(float[] array, float x){

        if (x<0) return array[0];
        else if (x>array.length-1) return array[array.length-1];
        else {
            int x0 = (int) Math.floor(x);
            return ((1-x+x0)*array[x0]+(x-x0)*array[x0+1]);}
    }

    // TODO: could be moved to NanoJ-Core
    public static double linearInterpolation(double[] array, float x){

        if (x<0) return array[0];
        else if (x>array.length-1) return array[array.length-1];
        else {
            int x0 = (int) Math.floor(x);
            return ((1-x+x0)*array[x0]+(x-x0)*array[x0+1]);}
    }

    // Adapted from nanoj.core.java.gui.ApplyDriftCorrection;
    // TODO: could be moved to NanoJ-Core
    class ThreadedRotateTranslate extends Thread {
        private final double angle;
        private final double shiftX;
        private final double shiftY;

        private final int slice;
        private final ImageStack ims;

        public ThreadedRotateTranslate(ImageStack ims, int slice, double shiftX, double shiftY, double angle) {
            this.shiftX = shiftX;
            this.shiftY = shiftY;
            this.angle = angle;
            this.ims = ims;
            this.slice = slice;
        }

        @Override
        public void run() {
            ImageProcessor ip = ims.getProcessor(slice).duplicate(); // TODO duplicate or not duplicate, that is the question????
            ip.setInterpolationMethod(ip.BICUBIC);
            ip.rotate(angle);
            ip.translate(shiftX, shiftY);
            ims.setProcessor(ip, slice);
        }
    }

    // Adapted from nanoj.core.java.gui.ApplyDriftCorrection;
    // TODO: could be moved to NanoJ-Core
    class ThreadedRotate extends Thread {
        private final double angle;
        private final int slice;
        private final ImageStack ims;

        public ThreadedRotate(ImageStack ims, int slice, double angle) {
            this.angle = angle;
            this.ims = ims;
            this.slice = slice;
        }

        @Override
        public void run() {
            ImageProcessor ip = ims.getProcessor(slice).duplicate(); // needs to be duplicated for rotate to work
            ip.setInterpolationMethod(ip.BICUBIC);
            ip.rotate(angle);
            ims.setProcessor(ip, slice);
        }
    }


}