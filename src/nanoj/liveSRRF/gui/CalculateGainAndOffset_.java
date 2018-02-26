package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.pow;

public class CalculateGainAndOffset_  implements PlugIn {

    @Override
    public void run(String s) {

        String dirPath = IJ.getDirectory("Choose directory to open...");
        if (dirPath == null) return;

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Simple Radiality");
        gd.addNumericField("Frames per calculation", 500, 0);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int nFrames = (int) gd.getNextNumber();

        //IJ.run("Image Sequence...", "open=["+dirPath+"] sort use");
        IJ.run("Image Sequence...", "open=["+dirPath+"] sort");
        ImagePlus imp = IJ.getImage();
        ImageStack ims = imp.getImageStack();

        assert(ims.getSize() % nFrames == 0);

        ImageStack[] imsCalculated = this.calculateMeanAndVarianceProjection(ims, nFrames);
        ImageStack imsMean = imsCalculated[0];
        ImageStack imsVar  = imsCalculated[1];

        // mean subtract and var subtract
        float[] mean0 = ((float[]) imsMean.getProcessor(1).getPixels()).clone();
        float[] var0  = ((float[]) imsVar.getProcessor(1).getPixels()).clone();

        for (int f=1; f<=imsMean.getSize(); f++) {
            IJ.showProgress(f, ims.getSize());
            IJ.showStatus("Mean and Var subtracting frame: "+f+"/"+ims.getSize());

            float[] mean1 = (float[]) imsMean.getProcessor(f).getPixels();
            float[] var1  = (float[]) imsVar.getProcessor(f).getPixels();

            for (int n=0; n<mean1.length; n++) {
                mean1[n] -= mean0[n];
                var1[n]  -= var0[n];
            }
        }

        new ImagePlus("Mean", imsMean).show();
        new ImagePlus("Var", imsVar).show();
    }

    private ImageStack[] calculateMeanAndVarianceProjection(ImageStack ims, int nFrames) {
        int w = ims.getWidth();
        int h = ims.getHeight();

        ImageStack imsMean = new ImageStack(w, h);
        ImageStack imsVar = new ImageStack(w, h);

        float[] mean = new float[w*h];
        float[] var  = new float[w*h];

        int counter = 1;

        for (int s=1; s<=ims.getSize(); s++) {
            IJ.showProgress(s, ims.getSize());
            IJ.showStatus("Processing frame: "+s+"/"+ims.getSize());

            FloatProcessor frame = ims.getProcessor(s).convertToFloatProcessor();
            float[] pixels = (float[]) frame.getPixels();

            for (int n=0; n<pixels.length; n++) {
                mean[n] += (pixels[n]-mean[n]) / counter;
                var[n] += pow(pixels[n]-mean[n], 2) / counter; // check on this
            }

            //counter ++;

            if (counter == nFrames) {
                FloatProcessor fpMean = new FloatProcessor(w, h, mean);
                FloatProcessor fpVar = new FloatProcessor(w, h, var);

                imsMean.addSlice(fpMean);
                imsVar.addSlice(fpVar);

                mean = new float[w*h];
                var  = new float[w*h];
                counter = 1;
            }
            else {
                counter++;
            }
        }

        return new ImageStack[] {imsMean, imsVar};
    }
}
