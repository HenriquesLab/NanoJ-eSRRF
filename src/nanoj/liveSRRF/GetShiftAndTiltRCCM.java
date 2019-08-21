package nanoj.liveSRRF;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core.java.gui.ApplyDriftCorrection;
import nanoj.core.java.image.transform.TranslateOrRotateImage;
import nanoj.core.java.tools.Log;
import static java.lang.Math.*;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.array.ArrayInitialization.initializeFloatAndGrowthFill;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getMaxFindByOptimization;

public class GetShiftAndTiltRCCM {
    private static final TranslateOrRotateImage translateOrRotateImage = new TranslateOrRotateImage();

    private ImageProcessor ipRef;
    private float[] angleArray;
    private float angleStep;
    private static final Log log = new Log();
    //    private int radiusX, radiusY; // TODO: include radius X and Y to speed things up?
    public ImageStack[] imsRCCMap;

    // Initialisation
    public GetShiftAndTiltRCCM(ImageProcessor ipRef, int nAngles, float maxAngle){
        this.ipRef = ipRef.duplicate();
        this.angleStep = 2*(float) toRadians(maxAngle)/(nAngles-1);
        this.angleArray = initializeFloatAndGrowthFill(nAngles, -(float) toRadians(maxAngle), angleStep);
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
            imsRCCMap[s - 1] = calculateCrossCorrelationMap(ipRef, imsRotated, true);
        }
    }

    // This method calculates the shift and tilt from the RCCM ImageStack array ----------------------------------------
    public double[][] getShiftAndTiltfromRCCM(){

        int width = imsRCCMap[0].getWidth();
        int nSlices = imsRCCMap.length;

        double[] shiftX = new double[nSlices];
        double[] shiftY = new double[nSlices];
        double[] theta = new double[nSlices];

        for(int s=0; s<nSlices; s++) {
            ImageStack imsRCCMapSlice = imsRCCMap[s];

            // calculate position and value of maximum pixel TODO: avoid spurious pixels that would skew this
            float v, vMax = 0;
            int xMax = 0, yMax = 0, aMax = 0;
            for (int a = 0; a<imsRCCMapSlice.getSize(); a++) {
                float[] pixelsSlice = (float[]) imsRCCMapSlice.getPixels(a+1);
                for (int p = 0; p < pixelsSlice.length; p++) {
                    v = pixelsSlice[p];
                    if (v > vMax) {
                        vMax = v;
                        xMax = p % width;
                        yMax = p / width;
                        aMax = a;
                    }
                }
            }
//            IJ.log("Slice: "+s);
//            IJ.log("aMax: "+aMax);
            theta[s] = toDegrees(angleArray[aMax]);
//            IJ.log("theta (in degrees): "+theta[s]);

            FloatProcessor fpRCCMslice = imsRCCMapSlice.getProcessor(aMax+1).convertToFloatProcessor();
            float[] shiftAndTilt = getMaxFindByOptimization(fpRCCMslice);

            shiftX[s] = (imsRCCMapSlice.getWidth()-1)/2 - shiftAndTilt[0];
            shiftY[s] = (imsRCCMapSlice.getHeight()-1)/2 - shiftAndTilt[1];
        }

        return new double[][] {shiftX, shiftY, theta};
    }

    public static ImagePlus applyCorrection(ImagePlus imp, double[] shiftX, double[] shiftY, double[] theta){

        float[] floatTheta = new float[theta.length];
        for (int i = 0 ; i < theta.length; i++) {
            floatTheta[i] = (float) toRadians(theta[i]);
        }

        ApplyDriftCorrection adc = new ApplyDriftCorrection();
        ImageStack ims = imp.getImageStack();
        ImageStack imsRotated = translateOrRotateImage.rotate(ims, floatTheta);

        ImagePlus impRot = new ImagePlus(imp.getShortTitle(), imsRotated);
        ImagePlus impAvgCorrected = adc.applyDriftCorrection(impRot, shiftX, shiftY);

        return impAvgCorrected;
    }


}