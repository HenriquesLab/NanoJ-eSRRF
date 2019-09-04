// Adapted from https://github.com/cmbruns/fiji-plugins/blob/master/Shannon_Entropy/Shannon_Entropy.java

package nanoj.liveSRRF;

import ij.*;
import ij.process.*;
import java.util.HashMap;

/** Compute the Shannon entropy of the current image.
 *  This statistical measure can be interpreted as the intrinsic bit depth of the image.
 *  Result in bits-per-pixel is shown in the log window.
 */

public class ShannonEntropyMeasure {

    public ImageStack ims;
    public double[] EntropyArray;

    public ShannonEntropyMeasure(ImageStack ims){
        this.ims = ims;
        int nSlices = ims.getSize();
        this.EntropyArray = new double[nSlices];


        HashMap<Integer, Integer> histogram = new HashMap<Integer, Integer>();

        // Count how many pixels there are of each intensity
        for (int n = 1; n <= nSlices; n++) {
            int totalPixels = 0;
            ImageProcessor ip = ims.getProcessor(n);
            for (int x = 0; x < ip.getWidth(); ++x) {
                for (int y = 0; y < ip.getHeight(); ++y) {
                    int val = ip.getPixel(x, y);
                    totalPixels++;
                    if (histogram.containsKey(val)) {
                        histogram.put(val, histogram.get(val) + 1);
                    } else {
                        histogram.put(val, 1);
                    }
                }
            }

            EntropyArray[n-1] = CalculateEntropy(histogram, totalPixels);
            IJ.showProgress(n, nSlices);
        }


    }


    public double CalculateEntropy(HashMap<Integer, Integer> histogram, int totalPixels){

        // Convert counts to probabilities and compute Shannon entropy
        double totalProbability = 0.0;
        double entropy = 0.0;
        int minIntensity = 1000000;
        int maxIntensity = -1;
        Integer modeIntensity = null;
        int modeCount = -1;

        for (Integer intensity : histogram.keySet()) {
            int count = histogram.get(intensity);
            // While we're here, store the mode statistic too
            if (count > modeCount) {
                modeCount = count;
                modeIntensity = intensity;
            }
            // And min and max
            if (minIntensity > intensity) {
                minIntensity = intensity;
            }
            if (maxIntensity < intensity) {
                maxIntensity = intensity;
            }
            double p = (double)count / (double)totalPixels;
            totalProbability += p;
            entropy -= p * Math.log(p) / Math.log(2.0); // bits
        }
        return entropy;
    }

}