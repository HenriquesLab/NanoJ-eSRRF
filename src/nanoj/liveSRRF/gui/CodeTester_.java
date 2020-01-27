package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;


public class CodeTester_ implements PlugIn {



    public void run(String arg) {

        float[][] pixelsMatrix = new float[3][10];
        pixelsMatrix[1][5] = 10;

        float[] pixels = pixelsMatrix[1];

        IJ.log("\\Clear");  // Clear the log window
        for (int i = 0; i < 10; i++) {
            IJ.log("Kabooom !! "+pixels[i]);
        }

    }
}
