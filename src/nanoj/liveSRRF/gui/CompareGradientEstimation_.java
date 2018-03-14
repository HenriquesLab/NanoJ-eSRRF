package nanoj.liveSRRF.gui;

import ij.*;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.CompareGradientEstimationCL;

public class CompareGradientEstimation_ implements PlugIn {

    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());

    @Override
    public void run(String arg) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        // initilizaing string for gradient estimation choice
        String[] GradMethods = new String[3];
        GradMethods[0] = "3-point gradient (classic)";
        GradMethods[1] = "Robert's cross local gradient";
        GradMethods[2] = "2-point local + interpolation";


        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Radial Gradient Convergence");
        gd.addChoice("Gradient estimation method", GradMethods, GradMethods[0]);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        String GradChosenMethod = gd.getNextChoice();

        IJ.log("Gradient method chosen: "+GradChosenMethod);

        ImageStack ims = imp.getImageStack();
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        CompareGradientEstimationCL gCL = new CompareGradientEstimationCL(w, h, GradChosenMethod);

        if (GradChosenMethod.equals("2-point local + interpolation")){
            w *= 2;
            h *= 2; }

        ImageStack imsGx = new ImageStack(w, h);
        ImageStack imsGy = new ImageStack(w, h);

        for (int n=1; n<=nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fpGxy = gCL.calculateGxGy(ims.getProcessor(n), 0, 0, GradChosenMethod);

//            if (GradChosenMethod.equals("2-point local + interpolation")){
//                imsGx.addSlice(InterpolateImage(fpGxy[0],2));
//                imsGy.addSlice(InterpolateImage(fpGxy[1],2));
//            }
//            else {
                imsGx.addSlice(fpGxy[0]);
                imsGy.addSlice(fpGxy[1]);
//            }



            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }

        new ImagePlus(imp.getTitle()+" - Gradient Gx", imsGx).show();
        new ImagePlus(imp.getTitle()+" - Gradient Gy", imsGy).show();


        //gCL.showPlatforms();
        gCL.release();
    }

    // Image interpolation function ------------------------------------------------------------------------------------
    public FloatProcessor InterpolateImage(FloatProcessor fp, int magnification){
        ImageProcessor ip = (ImageProcessor) fp;
        int width = ip.getWidth();
        int height = ip.getHeight();

        float[] PixelsInt = new float[width*magnification*height*magnification];

        for (int j=0; j<=height*magnification-1; j++) {
            for (int i=0; i<=width*magnification-1; i++) {
                int offset = i + j*width*magnification;
                double x = (double) i/magnification;
                double y = (double) j/magnification;
                PixelsInt[offset] = (float) ip.getBicubicInterpolatedPixel( x,  y, ip);
            }
        }

        FloatProcessor fpInt = new FloatProcessor(width*magnification, height*magnification, PixelsInt);
        return fpInt;
    }


}


