package nanoj.liveSRRF;

import ij.ImageStack;
import ij.process.FloatProcessor;

import static nanoj.core2.NanoJFHT.forwardFHT;
import static nanoj.core2.NanoJFHT.inverseFHT;

public class FourierInterpolation {

    public static FloatProcessor fhtSpaceInterpolation(FloatProcessor fp, int intFactor, boolean doMirrorPadding) {

        FloatProcessor fpFHT;
        if (doMirrorPadding){
            // Mirror padding of the original image
            FloatProcessor fpPadded = mirrorPadding(fp);
            fpFHT = forwardFHT(fpPadded);
        }
        else{
            fpFHT = forwardFHT(fp);
        }

        int w = fpFHT.getWidth();
        int h = fpFHT.getHeight();
        int wInt = w*intFactor;
        int hInt = h*intFactor;
        float[] pixelsFHT = (float[]) fpFHT.duplicate().getPixels();
        float[] pixelsFHTint = new float[wInt * hInt];
        float intF2 = (float) intFactor*intFactor;

        for (int p = 0; p < wInt; p++) {
            for (int q = 0; q < hInt; q++) {
                if (p <= (w / 2)) {
                    if (q <= (h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[p * w + q];
                    if ((q > (h / 2)) && (q <= (hInt - h / 2))) pixelsFHTint[p * wInt + q] = 0;
                    if (q > (hInt - h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[p*w + q + h - hInt];
                }

                if ((p > (w / 2)) && (p <= wInt - w / 2) && (q <= (h / 2))) pixelsFHTint[p * wInt + q] = 0;

                if (p > (wInt - w / 2)) {
                    if (q <= (h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[(p + w - wInt) * w + q];
                    if ((q > (h / 2)) && (q <= (hInt - h / 2))) pixelsFHTint[p * wInt + q] = 0;
                    if (q > (hInt - h / 2))
                        pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[(p + w - wInt) * w + q + h - hInt];
                }
            }
        }

        FloatProcessor fpFHTint = new FloatProcessor(wInt, hInt, pixelsFHTint);

        // Calculate the iFHT
        FloatProcessor fpInt = inverseFHT(fpFHTint, true);
        if (doMirrorPadding){
            int xROI = (fpFHTint.getWidth() - intFactor * fp.getWidth()) / 2;
            int yROI = (fpFHTint.getHeight() - intFactor * fp.getHeight()) / 2;
            fpInt.setRoi(xROI, yROI, intFactor * fp.getWidth(), intFactor * fp.getHeight());
            return fpInt.crop().convertToFloatProcessor();
        }
        else return fpInt;

    }

    public static ImageStack fhtSpaceInterpolation(ImageStack ims, int intFactor, boolean doMirrorPadding) {

        int nFrames = ims.getSize();
        ImageStack imsInt = new ImageStack();

        for (int i = 1; i <= nFrames; i++) {
            FloatProcessor fp = ims.getProcessor(i).duplicate().convertToFloatProcessor();
            FloatProcessor fpInt = fhtSpaceInterpolation(fp, intFactor, doMirrorPadding);
            imsInt.addSlice(fpInt);
        }

        return imsInt;

    }

    public static FloatProcessor mirrorPadding(FloatProcessor fp){

        int w = fp.getWidth();
        int h = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();
        float[] pixels3x3 = new float[3*w*3*h];

        int bitW = (int) Math.ceil(Math.log(w)/Math.log(2));
        if (Math.pow(2,bitW) == w) bitW++;
        bitW = (int) Math.pow(2, bitW);

        int bitH = (int) Math.ceil(Math.log(h)/Math.log(2));
        if (Math.pow(2,bitH) == h) bitH++;
        bitH = (int) Math.pow(2, bitH);

        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                pixels3x3[j*3*w + i] = pixels[(w-1-j)*w + (h-1-i)];
                pixels3x3[(j+h)*3*w + i] = pixels[j*w + (h-1-i)];
                pixels3x3[(j+2*h)*3*w + i] = pixels[(w-1-j)*w + (h-1-i)];

                pixels3x3[j*3*w + i+h] = pixels[(w-1-j)*w + i];
                pixels3x3[(j+h)*3*w + i+h] = pixels[j*w + i];
                pixels3x3[(j+2*h)*3*w + i+h] = pixels[(w-1-j)*w + i];

                pixels3x3[j*3*w + i+2*h] = pixels[(w-1-j)*w + (h-1-i)];
                pixels3x3[(j+h)*3*w + i+2*h] = pixels[j*w + (h-1-i)];
                pixels3x3[(j+2*h)*3*w + i+2*h] = pixels[(w-1-j)*w + (h-1-i)];
            }
        }

        FloatProcessor fpPadded = new FloatProcessor(3*w, 3*h, pixels3x3);

        int xROI = (3*w-bitW)/2; // TODO: currently works only if original image is even sized, check
        int yROI = (3*h-bitH)/2;

        fpPadded.setRoi(xROI, yROI, bitW, bitH);
        FloatProcessor fpCropped = fpPadded.crop().convertToFloatProcessor();

        return fpCropped;

    }

}
