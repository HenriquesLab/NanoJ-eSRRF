package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

import java.util.Random;

import static java.lang.Math.pow;

public class CalculateGainAndOffset_  implements PlugIn {

    public static boolean DEBUG = false;

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

        int w = ims.getWidth();
        int h = ims.getHeight();

        FloatProcessor fpGain = this.calculateGain(imsMean, imsVar);

        new ImagePlus("Gain", fpGain).show();
        new ImagePlus("Offset", new FloatProcessor(w,h,mean0)).show();

        if (DEBUG) {
            new ImagePlus("Mean", imsMean).show();
            new ImagePlus("Var", imsVar).show();
            showGainOffsetFit(imsMean, imsVar, fpGain, 10);
        }

    }

    // -----------------------------------------------------------------------------------------------------------------
    private ImageStack[] calculateMeanAndVarianceProjection(ImageStack ims, int nFrames) {
        int w = ims.getWidth();
        int h = ims.getHeight();

        ImageStack imsMean = new ImageStack(w, h);

        int counter = 1;
        float[] mean = new float[w*h];
        for (int s=1; s<=ims.getSize(); s++) {
            IJ.showProgress(s, ims.getSize());
            IJ.showStatus("Calculating mean for frame: "+s+"/"+ims.getSize());

            FloatProcessor frame = ims.getProcessor(s).convertToFloatProcessor();
            float[] pixels = (float[]) frame.getPixels();

            for (int n=0; n<pixels.length; n++) {
                mean[n] += pixels[n] / nFrames;
            }

            if (counter < nFrames) counter++;
            else {
                FloatProcessor fpMean = new FloatProcessor(w, h, mean);
                imsMean.addSlice(fpMean);
                mean = new float[w*h];
                counter = 1;
            }
        }
        //IJ.log("Mean frames="+imsMean.getSize());

        ImageStack imsVar = new ImageStack(w, h);

        counter = 1;
        float[] var  = new float[w*h];
        for (int s=1; s<=ims.getSize(); s++) {
            IJ.showProgress(s, ims.getSize());
            IJ.showStatus("Calculating var for frame: "+s+"/"+ims.getSize());

            FloatProcessor frame = ims.getProcessor(s).convertToFloatProcessor();
            FloatProcessor fpMean = (FloatProcessor) imsMean.getProcessor(imsVar.getSize()+1); // TODO: CHECK THIS IS CORRECT
            float[] pixels = (float[]) frame.getPixels();
            float[] pixelsMean = (float[]) fpMean.getPixels();

            for (int n=0; n<pixels.length; n++) {
                var[n] += pow(pixels[n] - pixelsMean[n], 2) / (nFrames - 1); // TODO: NEED TO CHECK BIAS VS UNBIAS
            }

            if (counter < nFrames) counter++;
            else {
                FloatProcessor fpVar = new FloatProcessor(w, h, var);
                imsVar.addSlice(fpVar);
                var = new float[w*h];
                counter = 1;
            }
        }

        return new ImageStack[] {imsMean, imsVar};
    }

    // -----------------------------------------------------------------------------------------------------------------
    private FloatProcessor calculateGain(ImageStack imsMean, ImageStack imsVar) { // TODO: have a toggle between the two fits with and without offset fit

        int w = imsMean.getWidth();
        int h = imsMean.getHeight();
        int nPixels = w*h;
        int nFrames = imsMean.getSize();

        float[] gain = new float[nPixels];
        double[] xy   = new double[nPixels];
        double[] xx   = new double[nPixels];

        for (int s=1; s<=nFrames; s++) {
            float[] pixelsMean = (float[]) imsMean.getProcessor(s).getPixels();
            float[] pixelsVar  = (float[]) imsVar.getProcessor(s).getPixels();

            for (int n=0; n<nPixels; n++) {
                xy[n] += pixelsMean[n]*pixelsVar[n];
                xx[n] += pixelsMean[n]*pixelsMean[n];
            }
        }

        for (int n=0; n<nPixels; n++) {
            gain[n] = (float) xy[n]/ (float) xx[n];
        }

        //ImageStack imsGain = new ImageStack(w, h);
        //imsGain.addSlice(new FloatProcessor(w, h, gain));

//        // first pass: compute means
//        float[] meanMean  = new float[nPixels];
//        float[] meanMean2 = new float[nPixels];
//        float[] meanVar   = new float[nPixels];
//
//        // first pass
//        for (int s=1; s<=nFrames; s++) {
//            float[] pixelsMean = (float[]) imsMean.getProcessor(s).getPixels();
//            float[] pixelsVar  = (float[]) imsVar.getProcessor(s).getPixels();
//
//            for (int n=0; n<nPixels; n++) {
//                meanMean[n] += pixelsMean[n] / nFrames;
//                meanMean2[n] += pixelsMean[n] * pixelsMean[n] / nFrames;
//                meanVar[n] += pixelsVar[n] / nFrames;
//            }
//        }
//
//        // x is mean and y is var
//        float[] xxbar = new float[nPixels];
//        float[] yybar = new float[nPixels];
//        float[] xybar = new float[nPixels];
//
//        // second pass
//        for (int s=1; s<=nFrames; s++) {
//            float[] pixelsMean = (float[]) imsMean.getProcessor(s).getPixels();
//            float[] pixelsVar  = (float[]) imsVar.getProcessor(s).getPixels();
//
//            for (int n=0; n<nPixels; n++) {
//                double deltaX = pixelsMean[n] - meanMean[n];
//                double deltaY = pixelsVar[n]  - meanVar[n];
//                xxbar[n] += deltaX * deltaX;
//                yybar[n] += deltaY * deltaY;
//                xybar[n] += deltaX * deltaY;
//            }
//        }
//
//        float[] gain = new float[nPixels];
//        float[] offset = new float[nPixels];
//
//        for (int n=0; n<nPixels; n++) {
//            gain[n] = xybar[n] / xxbar[n];
//            offset[n] = meanVar[n] - gain[n] * meanMean[n];
//        }
//
//        ImageStack imsGainOffset = new ImageStack(w, h);
//        imsGainOffset.addSlice(new FloatProcessor(w, h, gain));
//        imsGainOffset.addSlice(new FloatProcessor(w, h, offset));

        return new FloatProcessor(w, h, gain);

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

    // -----------------------------------------------------------------------------------------------------------------
    private void showGainOffsetFit(ImageStack imsMean, ImageStack imsVar, FloatProcessor fpGain, int nPlots) {
        int w = imsMean.getWidth();
        int h = imsMean.getHeight();
        int ns = imsMean.getSize();

        Random r = new Random();

        for (int i=0; i<nPlots; i++) {
            int x = r.nextInt(w);
            int y = r.nextInt(h);

            double[] xs = new double[ns];
            double[] ys0 = new double[ns];
            double[] ys1 = new double[ns];
            double[] gain = (double[]) fpGain.getPixels();


            for (int n=0; n<ns; n++) {
                xs[n] = (float) imsMean.getVoxel(x, y, n);
                ys0[n] = (float) imsVar.getVoxel(x, y, n);
                ys1[n] = gain[n]*xs[n];
            }
            Plot p = new Plot("Var vs Mean X="+x+" Y="+y, "Mean", "Var");
            p.addPoints(xs,ys1,Plot.LINE);
            p.addPoints(xs, ys0, Plot.CROSS);
            p.show();
        }
    }



}

