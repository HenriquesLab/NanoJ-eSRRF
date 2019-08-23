package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class JavaTest_ implements PlugIn {
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();
//
        ImageStack imsRot = imp.getStack().duplicate();

        for (int i = 0; i < imsRot.getSize(); i++) {
            ImageProcessor ip = imsRot.getProcessor(i+1);
            ip.setInterpolationMethod(ip.BICUBIC);
            ip.setBackgroundValue(0);
            ip.rotate(5.8d);
            imsRot.setProcessor(ip, i+1);
        }

        ImagePlus impRot = new ImagePlus("Rotated", imsRot);
        impRot.show();

    }
}
