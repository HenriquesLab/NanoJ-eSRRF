package nanoj.liveSRRF;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import static java.lang.Float.NaN;

public class OpticalFlowMagnitudeEstimator {

    float[] gradMag, pixelsRef;
    int w,h;
    float noiseStd;

    public OpticalFlowMagnitudeEstimator(ImageProcessor ipRef, float noiseStd){
        this.pixelsRef = (float[]) ipRef.convertToFloatProcessor().getPixels();
        this.gradMag = getGradientMag(ipRef);
        this.w = ipRef.getWidth();
        this.h = ipRef.getHeight();
        this.noiseStd = noiseStd;
    }

    public FloatProcessor calculateMagnitude(ImageProcessor ip){
        float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();
        float[] flowMag = new float[w*h];

        for (int i = 0; i < w*h; i++) {
            if (gradMag[i] > 2.5f*noiseStd) flowMag[i] = Math.abs(pixels[i] - pixelsRef[i])/gradMag[i];
            else flowMag[i] = NaN;
        }

        return new FloatProcessor(w,h,flowMag);
    }


    public static float[] getGradientMag(ImageProcessor ip){

        int w = ip.getWidth();
        int h = ip.getHeight();

        float[] gradMag = new float[w*h];
        float[] pixels = (float[]) ip.convertToFloatProcessor().getPixels();

        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                int x0 = Math.max(x-1,0);
                int x1 = Math.min(x+1, w-1);
                int y0 = Math.max(y-1,0);
                int y1 = Math.min(y+1, h-1);

                gradMag[x*w + y] = 0.5f * (float) Math.sqrt((pixels[x1*w+y]-pixels[x0*w+y])*(pixels[x1*w+y]-pixels[x0*w+y]) + (pixels[x*w+y1]-pixels[x*w+y0])*(pixels[x*w+y1]-pixels[x*w+y0]));
            }
        }

        return gradMag;
    }
}
