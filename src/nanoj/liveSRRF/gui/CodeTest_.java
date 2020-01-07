package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.plugin.PlugIn;

import ij.process.FloatProcessor;
import nanoj.liveSRRF.EntropyMeasure;


public class CodeTest_ implements PlugIn {
    public void run(String arg) {

        int imageSize = 64;
        float[] pixels = new float[imageSize*imageSize];

//        float period = 13.8f; // in pixels
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Generating sinewaves ~~~~~~~");
        gd.addNumericField("Period (in pixels): ", 13.8, 2);
        gd.showDialog();

        float period = (float) gd.getNextNumber();

        for (int y = 0; y < imageSize; y++) {
            for (int x = 0; x < imageSize; x++) {
                pixels[y*imageSize + x] = (float) Math.cos(2*Math.PI*x/period);
            }
        }

        ImageStack ims = new ImageStack(imageSize, imageSize);
        ims.addSlice(new FloatProcessor(imageSize, imageSize, pixels));
        ImagePlus imp = new ImagePlus("IMAGE", ims);
        imp.show();

    }


}
