package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import nanoj.liveSRRF.FHT3D_CL;

import static nanoj.liveSRRF.FHT3D_Utilities.*;


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
        gd.addCheckbox("Swap quadrants", false);
        gd.addCheckbox("Get FFT magnitude", false);
        gd.addCheckbox("Reshape for interpolation", false);
        gd.addNumericField("Magnification", 2, 0);
        gd.addCheckbox("Get FRC resolution", false);
        gd.addNumericField("Pixel size (nm): ", 100, 2);
        gd.showDialog();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        boolean inverse = gd.getNextBoolean();
        boolean swapQuadrants = gd.getNextBoolean();
        boolean getFFTmag = gd.getNextBoolean();
        boolean reshapeForInt = gd.getNextBoolean();
        int mag = (int) gd.getNextNumber();
        boolean getFRC = gd.getNextBoolean();
        float pixelSize = (float) gd.getNextNumber();

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

        if (swapQuadrants){
            time = System.currentTimeMillis();
            ImageStack imsFHTSwap = swapQuadrants3D(imsFHT);
            time = System.currentTimeMillis() - time;
            if (PROFILINGMSG) {
                IJ.log("Quadrant swapping took " + (time) + " ms.");
            }
            ImagePlus impFHTSwap = new ImagePlus("FHT stack - quadrant swapped", imsFHTSwap);
            impFHTSwap.show();
        }

        if (getFFTmag){
            time = System.currentTimeMillis();
            ImageStack imsFFTmag = getFFTmagnitude(imsFHT);
            time = System.currentTimeMillis() - time;
            if (PROFILINGMSG) {
                IJ.log("Calculating FFT mag. took " + (time) + " ms.");
            }
            ImagePlus impFFTmag = new ImagePlus("FHT stack - FFT magnitude", imsFFTmag);
            impFFTmag.show();
        }

        if (getFRC){
            ImagePlus imp2 = IJ.openImage();
            if (imp2 == null) return;
            imp2.show();

            f = new FHT3D_CL(w, h, d);
            time = System.currentTimeMillis();
            f.run(imp2.getImageStack(), inverse);
            time = System.currentTimeMillis() - time;
            if (PROFILINGMSG) {
                IJ.log("OpenCL execution took " + (time) + " ms.");
            }

            ImageStack imsFHT2 = FHT3D_CL.imsFHT;
            ImagePlus impFHT2 = new ImagePlus("FHT stack 2", imsFHT2);
            impFHT2.show();
            getFSCresolution(imsFHT, imsFHT2, pixelSize);
        }


        IJ.log("All done.");

    }

}
