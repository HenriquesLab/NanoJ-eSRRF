package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.liveSRRF.GaussianFitMinimizer;

public class EstimatePSFsize_ implements PlugIn {

    public void run(String arg) {


        IJ.log("\\Clear");  // Clear the log window
        IJ.log("Estimating PSF parameters");
        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Extract PSF parameters...");
        gd.addChoice("Method: ", new String[]{"Gaussian", "Integrated gaussian"}, "Integrated gaussian");
        gd.showDialog();

        String method = gd.getNextChoice();
        IJ.log("Chosen method: "+method);
        int model = 0; // default: Gaussian
        if (method.equals("Gaussian")) model = GaussianFitMinimizer.GAUSSIAN;
        if (method.equals("Integrated gaussian")) model = GaussianFitMinimizer.INTEGRATED_GAUSSIAN;


        IJ.log("Extracting X and Y integrated signal");
        FloatProcessor fp = imp.getProcessor().convertToFloatProcessor();
        double[] signalX = new double[fp.getWidth()];
        double[] signalY = new double[fp.getHeight()];

        for (int x = 0; x < fp.getWidth(); x++) {
            for (int y = 0; y < fp.getHeight(); y++) {
                signalX[x] += fp.getPixelValue(x,y);
                signalY[y] += fp.getPixelValue(x,y);
            }
        }

        IJ.log("Plotting the curves");
        Plot tracePlots = new Plot("X and Y profiles", "Pixel position", "Intensity (AU)");
        tracePlots.add("X profile",signalX);
        tracePlots.add("Y profile",signalY);
        tracePlots.show();

        double initialSigma = 1.0; // in pixels TODO: estimate initial Sigma better? 1 pixel seems sensible though.
        double initialX0 = (double) fp.getWidth()/2.0;
        double initialY0 = (double) fp.getWidth()/2.0;

        IJ.log("Getting minimizer ready...");
        GaussianFitMinimizer minimizer = new GaussianFitMinimizer(fp, initialSigma, initialX0, initialY0, model);
        double[] fitParameters = minimizer.calculate();


        IJ.log("X0: "+fitParameters[0]);
        IJ.log("Y0: "+fitParameters[1]);
        IJ.log("Sigma: "+fitParameters[2]);
        IJ.log("Amplitude: "+fitParameters[3]);
        IJ.log("Background: "+fitParameters[4]);


    }
}
