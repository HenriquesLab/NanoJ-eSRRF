package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;


public class CodeTester_ implements PlugIn {



    public void run(String arg) {

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Fluorescence simulator");
        gd.addNumericField("Number of frames: ", 10, 0);
        gd.addNumericField("Exposure time (ms): ", 10.0, 1);
        gd.showDialog();

        int nFrames = (int) gd.getNextNumber();
        double exposure = (double) gd.getNextNumber();
        exposure /= 1000f;

        double simDuration = (double) nFrames * exposure; // in seconds
        IJ.log("Acquisition time: "+simDuration+"s ("+nFrames+" frames at "+(exposure*1000)+"ms exposure)");

    }
}
