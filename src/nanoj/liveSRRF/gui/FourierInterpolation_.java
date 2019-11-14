package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import static nanoj.core2.NanoJFHT.*;

public class FourierInterpolation_ implements PlugIn { // TODO: check what happens when using any sized image
    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        GenericDialog gd = new GenericDialog("FHT interpolator");
        gd.addNumericField("Interpolation factor", 4, 0);
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        int interpolationFactor = (int) gd.getNextNumber();
        ImageStack ims = imp.getStack();
        int nFrames = ims.getSize();

        ImageStack imsPadded = new ImageStack();
        ImageStack imsFHT = new ImageStack();
        ImageStack imsFHTint = new ImageStack();
        ImageStack imsInt = new ImageStack();

        for (int i = 1; i <= nFrames; i++) {
            // Mirror padding of the original image
            FloatProcessor fp = ims.getProcessor(i).duplicate().convertToFloatProcessor();
            FloatProcessor fpPadded = mirrorPadding(fp);
            imsPadded.addSlice(fpPadded);

            // Calculate the FHT
            FloatProcessor fpFHT = forwardFHT(fpPadded);
            imsFHT.addSlice(fpFHT);

            // Calculate the corresponding FHT from interpolation
            FloatProcessor fpFHTint = fhtSpaceInterpolation(fpFHT, interpolationFactor);
            imsFHTint.addSlice(fpFHTint);

            // Calculate the iFHT
            FloatProcessor fpInt = inverseFHT(fpFHTint, true);
            int xROI = (fpFHTint.getWidth() - interpolationFactor * fp.getWidth()) / 2;
            int yROI = (fpFHTint.getHeight() - interpolationFactor * fp.getHeight()) / 2;
            fpInt.setRoi(xROI, yROI, interpolationFactor * fp.getWidth(), interpolationFactor * fp.getHeight());

            FloatProcessor fpIntCropped = fpInt.crop().convertToFloatProcessor();
            imsInt.addSlice(fpIntCropped);
        }

        ImagePlus impPadded = new ImagePlus("Mirror-padded image", imsPadded);
        ImagePlus impFHT = new ImagePlus("FHT", imsFHT);
        ImagePlus impFHTint = new ImagePlus("FHT-interpolated image (in FHT space)", imsFHTint);
        ImagePlus impInt = new ImagePlus("FHT-interpolated image", imsInt);

        impPadded.show();
        IJ.run(impPadded, "Enhance Contrast", "saturated=0.35");
        impFHT.show();
        IJ.run(impFHT, "Enhance Contrast", "saturated=0.35");
        impFHTint.show();
        IJ.run(impFHTint, "Enhance Contrast", "saturated=0.35");
        impInt.show();
        IJ.run(impInt, "Enhance Contrast", "saturated=0.35");

    }

    public static FloatProcessor fhtSpaceInterpolation(FloatProcessor fpFHT, int intFactor){

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

        return new FloatProcessor(wInt, hInt, pixelsFHTint);
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
