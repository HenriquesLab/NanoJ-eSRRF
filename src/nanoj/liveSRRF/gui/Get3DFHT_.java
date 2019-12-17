package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.liveSRRF.FHT3D_CL;


public class Get3DFHT_ implements PlugIn {

    private static boolean DEBUG = true;
    private static boolean PROFILINGMSG = true;


    @Override
    public void run(String s) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("3D FHT (OpenCL LOCI)");
        gd.addCheckbox("Inverse FHT", false);
        gd.addCheckbox("Reshape for interpolation", false);
        gd.addNumericField("Magnification", 2, 0);
        gd.showDialog();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        boolean inverse = gd.getNextBoolean();
        boolean reshapeForInt = gd.getNextBoolean();
        int mag = (int) gd.getNextNumber();

        long time; //used to track duration of ex
        int d = imp.getStackSize();
        int h = imp.getHeight();
        int w = imp.getWidth();

        ImageStack ims = imp.getStack();

        if (DEBUG) {
            IJ.log("Starting openCL setup...");
        }
        time = System.currentTimeMillis();
        FHT3D_CL f = new FHT3D_CL(w, h, d);
        time = System.currentTimeMillis() - time;
        if (PROFILINGMSG) {
            IJ.log("OpenCL setup took " + (time) + " ms.");
        }

        if (DEBUG) {
            IJ.log("Running openCL...");
        }
        time = System.currentTimeMillis();
        f.run(ims, inverse);
        time = System.currentTimeMillis() - time;
        if (PROFILINGMSG) {
            IJ.log("OpenCL execution took " + (time) + " ms.");
        }

        ImageStack imsFHT = FHT3D_CL.imsFHT;
        ImagePlus impFHT = new ImagePlus("FHT stack", imsFHT);
        impFHT.show();

        if (reshapeForInt) {
            time = System.currentTimeMillis();
            ImageStack imsFHTforInt = reshape3DFHTforInterpolation(imsFHT, mag);
            time = System.currentTimeMillis() - time;
            if (PROFILINGMSG) {
                IJ.log("Reshaping took " + (time) + " ms.");
            }
            ImagePlus impFHTforInt = new ImagePlus("FHT stack - reshaped", imsFHTforInt);
            impFHTforInt.show();
        }

        IJ.log("All done.");

    }


    private ImageStack reshape3DFHTforInterpolation(ImageStack ims, int magnification){

        int nSlices = ims.getSize();
        ImageStack imsOut = new ImageStack(ims.getWidth()*magnification, ims.getHeight()*magnification);
        FloatProcessor fpInt;
        for (int z = 0; z < nSlices*magnification; z++) {
            if (z < nSlices/2){
//                IJ.log("z position: "+z+" ---> mode 1");
                FloatProcessor fp = ims.getProcessor(z+1).convertToFloatProcessor();
                fpInt = fhtSpaceInterpolation(fp, magnification);
            }
            else if (z < (nSlices*magnification - nSlices/2)){
//                IJ.log("z position: "+z+" ---> mode 2");
                float[] pixels = new float[ims.getWidth()*magnification * ims.getHeight()*magnification];
                fpInt = new FloatProcessor(ims.getWidth()*magnification, ims.getHeight()*magnification, pixels);
            }
//            else if (z >= (nSlices*magnification - nSlices/2)){
            else {
//                IJ.log("z position: "+z+" ---> mode 3");
                int zOffset = z - nSlices*(magnification - 1);
//                IJ.log("z offset: "+zOffset);
                FloatProcessor fp = ims.getProcessor(zOffset+1).convertToFloatProcessor();
                fpInt = fhtSpaceInterpolation(fp, magnification);
            }

            imsOut.addSlice(fpInt);
        }

        return imsOut;

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
