package nanoj.liveSRRF;

import ij.ImagePlus;
import ij.ImageStack;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJProfiler;

import static nanoj.core.java.image.drift.EstimateShiftAndTilt.MAX_FITTING;
import static nanoj.core.java.image.transform.CrossCorrelationMap.calculateCrossCorrelationMap;
import static nanoj.core.java.image.drift.EstimateShiftAndTilt.getShiftFromCrossCorrelationPeak;

public class XYShiftCalculator {

    private ImageStack ims, // TODO: does this duplicate the memory used then ??
                         imsSubStack;
//    private final int radiusCCM = 20; // TODO: add this parameter as advanced parameter in GUI? Currently calculates the peak on the whole CCM maps

    public float[] shiftX, shiftY;

    // -- Constructor --
    public XYShiftCalculator(ImagePlus imp){
        this.ims = imp.getImageStack(); // This is the entire raw data
    }

    // -- Calculate shift Array using Cross-correlation matrix --
    public void calculateShiftArray(int indexStart, int nFrameForSRRF) {

        ImageProcessor ipRef = ims.getProcessor(indexStart);

        shiftX = new float[nFrameForSRRF];
        shiftY = new float[nFrameForSRRF];

        imsSubStack = new ImageStack(ims.getWidth(), ims.getHeight());

        for (int s = 0; s < nFrameForSRRF; s++) {
            imsSubStack.addSlice(ims.getProcessor(indexStart+s));
        }

        ImageStack imsCMM = calculateCrossCorrelationMap(ipRef, imsSubStack, true);
        float[][] shiftArray = getShiftFromCrossCorrelationPeak(imsCMM, MAX_FITTING);

        for (int s = 0; s < nFrameForSRRF; s++) {
            shiftX[s] = -shiftArray[1][s];
            shiftY[s] = -shiftArray[2][s];
            System.out.format("Frame=%d shiftX=%f shiftY=%f\n",s,shiftX[s],shiftY[s]);
        }

    }

}
