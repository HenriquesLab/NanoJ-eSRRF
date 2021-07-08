package nanoj.liveSRRF;

import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import nanoj.core.java.array.ArrayMath;

public class EntropyMeasure {

    private final int nBins = 100000; // number of bins
    public double[] entropyArray;

    public EntropyMeasure(ImageStack ims){

        int nSlices = ims.getSize();
        entropyArray = new double[nSlices];
        float[] minMax = getMinAndMaxFromImageStack(ims);
        IJ.log("Min="+minMax[0]+" and max="+minMax[1]);

        for (int s = 0; s < nSlices; s++) {
            entropyArray[s] = 0;
            int[] counts = new int[nBins];
            FloatProcessor fp = ims.duplicate().getProcessor(s+1).convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();
//            IJ.log("This pixel: "+pixels[0]+" and this one: "+pixels[10]);
            // Build the histogram
            for (int i = 0; i < ims.getWidth()*ims.getHeight(); i++) {
                counts[(int) ((nBins-1)*(pixels[i]-minMax[0])/(minMax[1]-minMax[0]))]++;
            }
            int sumCounts = sumIntArray(counts);
//            IJ.log("Sum of counts: "+sumCounts);
            // Compute the information entropy
            for (int i = 0; i < ims.getWidth()*ims.getHeight(); i++) {
                if (counts[i] > 0) entropyArray[s] -= (double) counts[i]/(double) sumCounts * Math.log((double) counts[i]/(double) sumCounts)/Math.log(2.0d);
            }
        }

    }


    // Get min and max from the whole imageStack
    private float[] getMinAndMaxFromImageStack(ImageStack ims){
        float min = Float.MAX_VALUE;
        float max = -Float.MAX_VALUE;
        for (int s = 0; s < ims.getSize(); s++) {
            FloatProcessor fp = ims.duplicate().getProcessor(s+1).convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();
            for (int i = 0; i < ims.getWidth()*ims.getHeight(); i++) {
                if (pixels[i]>max) max = pixels[i];
                if (pixels[i]<min) min = pixels[i];
            }
        }

        return new float[]{min, max};
    }


    private int sumIntArray(int[] counts){
        int sum = 0;
        for (int i = 0; i < counts.length; i++) {
            sum += counts[i];
        }
        return sum;
    }

}
