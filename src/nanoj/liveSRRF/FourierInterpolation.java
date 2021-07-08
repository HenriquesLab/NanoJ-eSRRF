package nanoj.liveSRRF;

import ij.ImageStack;
import ij.process.FloatProcessor;

import static java.lang.Integer.max;
import static nanoj.core2.NanoJFHT.forwardFHT;
import static nanoj.core2.NanoJFHT.inverseFHT;

public class FourierInterpolation {

    public static FloatProcessor fhtSpaceInterpolation(FloatProcessor fp, int intFactor, boolean doMirrorPadding) {

        FloatProcessor fpFHT;
        if (doMirrorPadding){
            // Mirror padding of the original image
            FloatProcessor fpPadded = mirrorPaddingEvenSquare(fp);
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


        int pInt, qInt;
        for (int q = 0; q < h; q++) {
            for (int p = 0; p < w; p++) {
                if (p < w/2) pInt = p;
                else pInt = p + (intFactor-1)*w;
                if (q < h/2) qInt = q;
                else qInt = q + (intFactor-1)*h;

                pixelsFHTint[pInt + qInt*wInt] = pixelsFHT[p + q*w];
            }
        }


        FloatProcessor fpFHTint = new FloatProcessor(wInt, hInt, pixelsFHTint);
//        return fpFHTint;

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

//    public static FloatProcessor mirrorPadding(FloatProcessor fp){
//
//        int w = fp.getWidth();
//        int h = fp.getHeight();
//        float[] pixels = (float[]) fp.getPixels();
//        float[] pixels3x3 = new float[3*w*3*h];
//
//        for (int j = 0; j < h; j++) {
//            for (int i = 0; i < w; i++) {
//                pixels3x3[j*3*w + i] = pixels[(h-1-j)*w + (w-1-i)]; // Top left
//                pixels3x3[(j+h)*3*w + i] = pixels[j*w + (w-1-i)];  // Middle left
//                pixels3x3[(j+2*h)*3*w + i] = pixels[(h-1-j)*w + (w-1-i)]; //Bottom Left
//
//                pixels3x3[j*3*w + i+w] = pixels[(h-1-j)*w + i];
//                pixels3x3[(j+h)*3*w + i+w] = pixels[j*w + i];
//                pixels3x3[(j+2*h)*3*w + i+w] = pixels[(h-1-j)*w + i];
//
//                pixels3x3[j*3*w + i+2*w] = pixels[(h-1-j)*w + (w-1-i)];
//                pixels3x3[(j+h)*3*w + i+2*w] = pixels[j*w + (w-1-i)];
//                pixels3x3[(j+2*h)*3*w + i+2*w] = pixels[(h-1-j)*w + (w-1-i)];
//            }
//        }
//
//        FloatProcessor fpPadded = new FloatProcessor(3*w, 3*h, pixels3x3);
//
//        int bitW = (int) Math.ceil(Math.log(w)/Math.log(2));
//        if (Math.pow(2,bitW) == w) bitW++;
//        bitW = (int) Math.pow(2, bitW);
//
//        int bitH = (int) Math.ceil(Math.log(h)/Math.log(2));
//        if (Math.pow(2,bitH) == h) bitH++;
//        bitH = (int) Math.pow(2, bitH);
//
//        int xROI = (3*w-bitW)/2; // TODO: currently works only if original image is even sized, check
//        int yROI = (3*h-bitH)/2;
//
//        fpPadded.setRoi(xROI, yROI, bitW, bitH);
//        FloatProcessor fpCropped = fpPadded.crop().convertToFloatProcessor();
//
//        return fpCropped;
//    }


    public static FloatProcessor mirrorPaddingEvenSquare(FloatProcessor fp){

        int w = fp.getWidth();
        int h = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();

        // Get the minimum power of two image in each dimension
        int bitW = (int) Math.ceil(Math.log(w)/Math.log(2));
//        if (Math.pow(2,bitW) == w) bitW++;
        bitW = (int) Math.pow(2, bitW);
        int bitH = (int) Math.ceil(Math.log(h)/Math.log(2));
//        if (Math.pow(2,bitH) == h) bitH++;
        bitH = (int) Math.pow(2, bitH);

//        IJ.log("bitH: "+bitH);
//        IJ.log("bitW: "+bitW);

        // Choose the max to make it square
        int outputDim = max(bitH, bitW);

        int scaleW = (int) Math.ceil((double) outputDim/w);
        if (scaleW%2 == 0) scaleW += 1; // make sure it's odd
        int scaleH = (int) Math.ceil((double) outputDim/h);
        if (scaleH%2 == 0) scaleH += 1;

//        IJ.log("scaleH: "+scaleH);
//        IJ.log("scaleW: "+scaleW);

        float[] pixelsPadded = new float[scaleW*w*scaleH*h];
        int p,q;

        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {

                for (int sH = 0; sH < scaleH ; sH++) {
                    for (int sW = 0; sW < scaleW ; sW++) {
                        if ((sW-(scaleW-1)/2)%2 == 0) p = i;
                        else p = w-1-i; // flip horizontally
                        if ((sH-(scaleH-1)/2)%2 == 0) q = j;
                        else q = h-1-j; // flip vertically
                        pixelsPadded[(j + sH*h)*scaleW*w + i + sW*w] = pixels[q*w + p];
                    }
                }
            }
        }

        FloatProcessor fpPadded = new FloatProcessor(scaleW*w, scaleH*h, pixelsPadded);
        int xROI = (scaleW*w-outputDim)/2;
        int yROI = (scaleH*h-outputDim)/2;

        fpPadded.setRoi(xROI, yROI, outputDim, outputDim);
        FloatProcessor fpCropped = fpPadded.crop().convertToFloatProcessor();

        return fpCropped;
    }

}
