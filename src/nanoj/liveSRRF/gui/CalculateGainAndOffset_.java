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

<<<<<<< HEAD
            ImagePlus imp = IJ.openImage(fileName);

=======
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
>>>>>>> ad46e023ec0e40644b42f494717c0e3afe9176eb
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

    private void calculateGain(ImageStack imsMean, ImageStack imsVar) {

        int w = imsMean.getWidth();
        int h = imsMean.getHeight();
        int nPixels = w*h;
        int nFrames = imsMean.getSize();

        // first pass: read in data, compute xbar and ybar

        double[] meanMean  = new double[nPixels];
        double[] meanMean2 = new double[nPixels];
        double[] meanVar   = new double[nPixels];

        // first pass
        for (int s=1; s<=nFrames; s++) {
            float[] pixelsMean = (float[]) imsMean.getProcessor(s).getPixels();
            float[] pixelsVar  = (float[]) imsVar.getProcessor(s).getPixels();

            for (int n=0; n<nPixels; n++) {
                meanMean[n] += pixelsMean[n] / nFrames;
                meanMean2[n] += pixelsMean[n] * pixelsMean[n] / nFrames;
                meanVar[n] += pixelsVar[n] / nFrames;
            }
        }

        // x is mean and y is var
        double[] xxbar = new double[nPixels];
        double[] yybar = new double[nPixels];
        double[] xybar = new double[nPixels];

        // second pass
        for (int s=1; s<=nFrames; s++) {
            float[] pixelsMean = (float[]) imsMean.getProcessor(s).getPixels();
            float[] pixelsVar  = (float[]) imsVar.getProcessor(s).getPixels();

            for (int n=0; n<nPixels; n++) {
                double deltaX = pixelsMean[n] - meanMean[n];
                double deltaY = pixelsVar[n]  - meanVar[n];
                xxbar[n] += deltaX * deltaX;
                yybar[n] += deltaY * deltaY;
                xybar[n] += deltaX * deltaY;
            }
        }

        double[] gain = new double[nPixels];
        double[] offset = new double[nPixels];

        for (int n=0; n<nPixels; n++) {
            gain[n] = xybar[n] / xxbar[n];
            offset[n] = meanVar[n] - gain[n] * meanMean[n];
        }

//        int MAXN = 1000;
//        int n = 0;
//        double[] x = new double[MAXN];
//        double[] y = new double[MAXN];
//
//        // first pass: read in data, compute xbar and ybar
//        double sumx = 0.0, sumy = 0.0, sumx2 = 0.0;
//        while(!StdIn.isEmpty()) {
//            x[n] = StdIn.readDouble();
//            y[n] = StdIn.readDouble();
//            sumx  += x[n];
//            sumx2 += x[n] * x[n];
//            sumy  += y[n];
//            n++;
//        }
//        double xbar = sumx / n;
//        double ybar = sumy / n;
//
//        // second pass: compute summary statistics
//        double xxbar = 0.0, yybar = 0.0, xybar = 0.0;
//        for (int i = 0; i < n; i++) {
//            xxbar += (x[i] - xbar) * (x[i] - xbar);
//            yybar += (y[i] - ybar) * (y[i] - ybar);
//            xybar += (x[i] - xbar) * (y[i] - ybar);
//        }
//        double beta1 = xybar / xxbar;
//        double beta0 = ybar - beta1 * xbar;
//
//        // print results
//        StdOut.println("y   = " + beta1 + " * x + " + beta0);
//
//        // analyze results
//        int df = n - 2;
//        double rss = 0.0;      // residual sum of squares
//        double ssr = 0.0;      // regression sum of squares
//        for (int i = 0; i < n; i++) {
//            double fit = beta1*x[i] + beta0;
//            rss += (fit - y[i]) * (fit - y[i]);
//            ssr += (fit - ybar) * (fit - ybar);
//        }
//        double R2    = ssr / yybar;
//        double svar  = rss / df;
//        double svar1 = svar / xxbar;
//        double svar0 = svar/n + xbar*xbar*svar1;
//        StdOut.println("R^2                 = " + R2);
//        StdOut.println("std error of beta_1 = " + Math.sqrt(svar1));
//        StdOut.println("std error of beta_0 = " + Math.sqrt(svar0));
//        svar0 = svar * sumx2 / (n * xxbar);
//        StdOut.println("std error of beta_0 = " + Math.sqrt(svar0));
//
//        StdOut.println("SSTO = " + yybar);
//        StdOut.println("SSE  = " + rss);
//        StdOut.println("SSR  = " + ssr);

    }
}
