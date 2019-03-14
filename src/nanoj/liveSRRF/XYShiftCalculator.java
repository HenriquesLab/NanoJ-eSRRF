package nanoj.liveSRRF;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.PointRoi;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJProfiler;

import static java.lang.Math.min;
import static nanoj.core2.NanoJCrossCorrelation.calculateCrossCorrelationMap;

public class XYShiftCalculator {

    private ImagePlus impCCM = null;
    private ImageStack ims; // TODO: does this duplicate the memory used then ??
    private final int radiusCCM = 8; // TODO: add this parameter as advanced parameter?

    public float[] shiftX, shiftY;

    private NanoJProfiler prof;

    // -- Constructor --
    public XYShiftCalculator(ImagePlus imp, NanoJProfiler prof){
        this.ims = imp.getImageStack();
        this.prof = prof;
    }


    // -- Calculate shift Array using Cross-correlation matrix --
    public boolean calculateShiftArray(int indexStart, int nFrameForSRRF) {

        ImageProcessor ipRef = ims.getProcessor(indexStart);
        ImageProcessor ipData;

        shiftX = new float[nFrameForSRRF];
        shiftY = new float[nFrameForSRRF];

        for (int s = 0; s < nFrameForSRRF; s++) {

            // Check if user is cancelling calculation
            IJ.showProgress(s, nFrameForSRRF);
            if (IJ.escapePressed()) {
                IJ.resetEscape();
                IJ.log("-------------------------------------");
                IJ.log("Reconstruction aborted by user.");
                return true;
            }

            // Grab the new frame from the list
            ipData = ims.getProcessor(s+1);

            // Estimate vibrations
            int id = prof.startTimer();
            float[] shift = calculateShift(ipRef, ipData);
            shiftX[s] = shift[0];
            shiftY[s] = shift[1];

            System.out.println("Frame=" + s + " shiftX=" + shiftX[s] + " shiftY=" + shiftY[s]);
            prof.recordTime("Drift Estimation", prof.endTimer(id));
        }

        return false;

    }


    ///////////////////////////////////////////
    /// Functions and other helping buddies ///
    ///////////////////////////////////////////

    // --- Calculate shift using Cross-correlation matrix ---
    private float[] calculateShift(ImageProcessor ipRef, ImageProcessor ipData) {

        FloatProcessor fpCCM = (FloatProcessor) calculateCrossCorrelationMap(ipRef, ipData, false);

        int windowSize = radiusCCM * 2 + 1;
        int xStart = fpCCM.getWidth() / 2 - radiusCCM;
        int yStart = fpCCM.getHeight() / 2 - radiusCCM;
        fpCCM.setRoi(xStart, yStart, windowSize, windowSize);
        fpCCM = (FloatProcessor) fpCCM.crop();

        double vMax = -Double.MAX_VALUE;
        double vMin = Double.MAX_VALUE;
        double xMax = 0;
        double yMax = 0;

        // first do coarse search for max
        for (int y = 1; y < windowSize - 1; y++) {
            for (int x = 1; x < windowSize - 1; x++) {
                double v = fpCCM.getf(x, y);
                if (v > vMax) {
                    vMax = v;
                    xMax = x;
                    yMax = y;
                }
                vMin = min(v, vMin);
            }
        }
        //System.out.println("xMax="+xMax+" yMax="+yMax);

        //vMax = -Double.MAX_VALUE;
        // do fine search for max
        for (double y = yMax; y < yMax + 1; y += 0.01) {
            for (double x = xMax; x < xMax + 1; x += 0.01) {
                double v = fpCCM.getBicubicInterpolatedPixel(x, y, fpCCM);
                if (v > vMax) {
                    vMax = v;
                    xMax = x;
                    yMax = y;
                }
            }
        }

        // recenter pixels
        float shiftX = (float) xMax - radiusCCM;
        float shiftY = (float) yMax - radiusCCM;

//        if (impCCM == null) {
//            impCCM = new ImagePlus("CCM Vibration Stabilisation", fpCCM);
//            impCCM.show();
//        }
//
//        impCCM.setProcessor(fpCCM);
//        impCCM.setRoi(new PointRoi(xMax + .5, yMax + .5));
//        impCCM.setDisplayRange(vMin, vMax);

        return new float[]{shiftX, shiftY};
    }

}
