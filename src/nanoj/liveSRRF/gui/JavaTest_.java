package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import static nanoj.liveSRRF.gui.GetSpatialCalibrationMFMdata_.getSortedIndices;

public class JavaTest_ implements PlugIn {
    public void run(String s) {

//        ImagePlus imp = WindowManager.getCurrentImage();
//        if (imp == null) imp = IJ.openImage();
//        imp.show();
//
//        ImageStack ims = imp.getImageStack();
//        ImageProcessor ip = ims.getProcessor(1);
//        FloatProcessor fp = ip.convertToFloatProcessor();
//
//        fp.setInterpolationMethod(ImageProcessor.BICUBIC);
//        IJ.log("0/0 "+fp.getInterpolatedPixel(0, 0));
//        IJ.log("1/1 "+fp.getInterpolatedPixel(1, 1));
//        IJ.log("0.5/0.5 "+fp.getInterpolatedPixel(0.5, 0.5));

        double[] array = {2,3,4,-1,0,1,-4,-3,-2};
        int[] sortedIndices = getSortedIndices(array);

        for (int i=0; i<array.length; i++){
            IJ.log("Array: "+array[i]);
        }
        IJ.log("----------");

        for (int i=0; i<array.length; i++){
            IJ.log("Array: "+sortedIndices[i]);
        }


    }
}
