package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core.java.image.calculator.FloatProcessorCalculator;
import scala.Int;


public class CodeTester_ implements PlugIn {

    public void run(String arg) {

        FloatProcessorCalculator fpc = new FloatProcessorCalculator();
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        ImageStack ims = imp.getStack();
        int nFrames = ims.getSize();
        IJ.log("# frames: "+nFrames);

        ImageStack imsMean = new ImageStack(ims.getWidth(), ims.getHeight());

        FloatProcessor fpMean = ims.getProcessor(1).duplicate().convertToFloatProcessor();
        imsMean.addSlice(fpMean.duplicate().convertToFloatProcessor());

        for (int i = 1; i < nFrames; i++) {
            fpMean = fpc.add(fpMean.duplicate().convertToFloatProcessor(), ims.getProcessor(i+1).duplicate().convertToFloatProcessor());
            FloatProcessor fpTemp = fpMean.duplicate().convertToFloatProcessor();
            fpTemp.multiply(1.0d/(double)(i+1));
            imsMean.addSlice(fpTemp);
        }

        fpMean.multiply(1.0d/(double)(nFrames));
        ImagePlus impMean = new ImagePlus("Mean image", imsMean);
        impMean.show();

        FloatProcessor fpSOFItau1 = new FloatProcessor(ims.getWidth(), ims.getHeight());
        ImageStack imsSOFItau1 = new ImageStack(ims.getWidth(), ims.getHeight());
        for (int i = 0; i < (nFrames-1); i++) {
            FloatProcessor fp1 = fpc.subtract(ims.getProcessor(i+1).duplicate().convertToFloatProcessor(), fpMean);
            FloatProcessor fp2 = fpc.subtract(ims.getProcessor(i+2).duplicate().convertToFloatProcessor(), fpMean);

            fpSOFItau1 = fpc.add(fpSOFItau1.duplicate().convertToFloatProcessor(), fpc.multiply(fp1, fp2));

            FloatProcessor fpTemp = fpSOFItau1.duplicate().convertToFloatProcessor();
            fpTemp.multiply(1.0d/(double)(i+1));
            imsSOFItau1.addSlice(fpTemp);
        }

        ImagePlus impSOFItau1 = new ImagePlus("SOFI (dTau = 1) cumulative stack", imsSOFItau1);
        impSOFItau1.show();



        FloatProcessor fpSOFItau0 = new FloatProcessor(ims.getWidth(), ims.getHeight());
        ImageStack imsSOFItau0 = new ImageStack(ims.getWidth(), ims.getHeight());
        for (int i = 0; i < (nFrames-1); i++) {
            FloatProcessor fp1 = fpc.subtract(ims.getProcessor(i+1).duplicate().convertToFloatProcessor(), fpMean);
            FloatProcessor fp2 = fpc.subtract(ims.getProcessor(i+1).duplicate().convertToFloatProcessor(), fpMean);

            fpSOFItau0 = fpc.add(fpSOFItau0.duplicate().convertToFloatProcessor(), fpc.multiply(fp1, fp2));

            FloatProcessor fpTemp = fpSOFItau0.duplicate().convertToFloatProcessor();
            fpTemp.multiply(1.0d/(double)(i+1));
            imsSOFItau0.addSlice(fpTemp);
        }

        ImagePlus impSOFItau0 = new ImagePlus("SOFI (dTau = 0) cumulative stack", imsSOFItau0);
        impSOFItau0.show();




    }
}
