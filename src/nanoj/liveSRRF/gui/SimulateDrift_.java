package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;

import java.awt.*;

public class SimulateDrift_ implements PlugIn {

    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();


        int nSlices = imp.getImageStack().getSize();
        int width = imp.getImageStack().getWidth();
        int height = imp.getImageStack().getHeight();
        String imageTitle = imp.getTitle();

//        String[] curveType = new String[2];
//        curveType[0] = "linear";
//        curveType[1] = "quadratic";


        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Drift simulation");
        gd.addMessage("-=-= x Drift simulation =-=-\n", headerFont);
//        gd.addChoice("Drift curve ", curveType, curveType[0]);
        gd.addNumericField("a1 (default: 0.1)", -0.2, 3);
        gd.addNumericField("a2 (default: 0.001)", 0.00,3);

        gd.addMessage("-=-= y Drift simulation =-=-\n", headerFont);
//        gd.addChoice("Drift curve ", curveType, curveType[0]);
        gd.addNumericField("a1 (default: 0.1)", -0.1, 3);
        gd.addNumericField("a2 (default: 0.001)", 0.001,3);

        gd.addHelp("https://www.youtube.com/watch?v=iuJDhFRDx9M");
        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            return;
        }

//        String curveTypeX = gd.getNextChoice();
        float a1X = (float) gd.getNextNumber();
        float a2X = (float) gd.getNextNumber();

//        String curveTypeY = gd.getNextChoice();
        float a1Y = (float) gd.getNextNumber();
        float a2Y = (float) gd.getNextNumber();

        ImagePlus impOutput = imp.duplicate();
        impOutput.setTitle(imageTitle + " - Drifted");

        float[] x = new float[nSlices];
        float[] y = new float[nSlices];
        float[] t = new float[nSlices];

        float sf;
        for (int s = 1; s < nSlices; s++) {
            sf = (float) s;
            t[s] = s;
            x[s] = a2X * sf * sf + a1X * sf + (float) Math.random();
            y[s] = a2Y * sf * sf + a1Y * sf + (float) Math.random();
        }

        Plot plot = new Plot("x/y drift", "Time (frame)", "Shift (pixel)");
        plot.setColor(Color.black);
        plot.addPoints(t, x, Plot.LINE);
        plot.setColor(Color.red);
        plot.addPoints(t, y, Plot.LINE);

        plot.addLegend("x drift\ny drift\n", "Top-Left");
        plot.setLimits(0.0,nSlices,-10,10);
        plot.show();


        for (int s = 1; s < nSlices; s++) {
            impOutput.setSlice(s+1);
            IJ.run(impOutput, "Translate...", "x="+x[s]+" y="+y[s]+" interpolation=Bicubic slice");
        }

        impOutput.show();

    }
}
