package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.plugin.Binner;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


public class CodeTester_ implements PlugIn {

    private ImagePlus imp;


    public void run(String arg) {


        // Get raw data
        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        FloatProcessor fp = imp.getProcessor().convertToFloatProcessor();
        FloatProcessor fpResized = fp.resize(fp.getWidth()/2, fp.getHeight()/2).convertToFloatProcessor();

        ImagePlus impResized = new ImagePlus("Image resized",fpResized);
        impResized.show();

        FloatProcessor fpScaled = fp.duplicate().convertToFloatProcessor();
        fpScaled.setInterpolationMethod(ImageProcessor.NONE);
        fpScaled.scale(0.5, 0.5);

        ImagePlus impScaled = new ImagePlus("Image scaled",fpScaled);
        impScaled.show();

        ImageProcessor ipBinned = fp.duplicate();
        Binner binner = new Binner();
        ipBinned = binner.shrink(ipBinned, 2,2, Binner.SUM);

        ImagePlus impBinned = new ImagePlus("Image binned", ipBinned);
        impBinned.show();


    }
}
