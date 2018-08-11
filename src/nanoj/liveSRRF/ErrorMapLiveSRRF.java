package nanoj.liveSRRF;


import ij.IJ;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.univariate.*;

import java.text.DecimalFormat;
import java.util.ArrayList;

import static java.lang.Math.*;
import static nanoj.core.java.array.ArrayCasting.toArray;
import static nanoj.core.java.array.ArrayMath.calculateMSE;
import static nanoj.core.java.array.ArrayMath.calculatePPMCC;

public class ErrorMapLiveSRRF {

    // Basic formats
    private static final String OVERBLUR_MESSAGE = "RSF constrained as no good minimum found";
    private static final float ROOT2 = (float) Math.sqrt(2);
    int width,
            height;

    float maxSigmaBoundary,
            alpha,
            beta,
            sigma_linear;

    double globalRMSE,
            globalPPMCC;

    float[] pixelsRef, ones;

    // Image formats
    FloatProcessor fpSR;

    // Other formats
    DecimalFormat df = new DecimalFormat("00.00");


    // -- Do something useful ---
    public void optimise(ImageProcessor impRef, ImageProcessor impSR, int magnification) {

        this.width = impRef.getWidth();
        this.height = impRef.getHeight();
        int nPixels = width * height;

        FloatProcessor fpRef = impRef.convertToFloatProcessor();
        // Set up pixel arrays for optimization
        this.pixelsRef = (float[]) fpRef.getPixels();
        this.ones = new float[nPixels];
        for(int i=0; i<nPixels; i++){ones[i] = 1;}

        this.maxSigmaBoundary = (4 / 2.35482f) * magnification; // this assumes Nyquist sampling in the ref image

        boolean overblurFlag = false;

        // Get SR FloatProcessor
        this.fpSR = impSR.convertToFloatProcessor();

        // UNIVARIATE OPTIMIZER - LINEAR MATCHING
        /// setup optimizer
        sigmaOptimiseFunction f = new sigmaOptimiseFunction(fpSR, pixelsRef, ones);
        UnivariateOptimizer optimizer = new BrentOptimizer(1e-10, 1e-14);
        /// run optimizer
        UnivariatePointValuePair result = optimizer.optimize(new MaxEval(1000),
                new UnivariateObjectiveFunction(f), GoalType.MINIMIZE, new SearchInterval(0, maxSigmaBoundary)); //NYQUIST ASSUMED

        this.sigma_linear = (float) result.getPoint();
        IJ.log("Best sigma is: " + sigma_linear);
        IJ.log("Best error is: " + result.getValue());

//        float[] errorList = toArray(f.getErrorList(), 1.0f);
//        float[] sigmaList = toArray(f.getSigmaList(), 1.0f);

//        Plot plot = new Plot("Brent optimiser: Error vs Sigma", "Sigma", "Error");
//        plot.addPoints(sigmaList, errorList, Plot.CIRCLE);
//        plot.show();

        if (abs(sigma_linear - maxSigmaBoundary) < 0.0001f) {
            overblurFlag = true;
            IJ.log(OVERBLUR_MESSAGE);
        }

        // GET ALPHA AND BETA
        FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
        FloatProcessor blurredOnes = new FloatProcessor(width, height, ones);
        blurredFp.blurGaussian(sigma_linear);
        blurredOnes.blurGaussian(sigma_linear);

        float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());
        this.alpha = aB[0];
        this.beta = aB[1];

        IJ.log("Alpha is: " + alpha + ", beta is: " + beta);
    }


    public FloatProcessor calculateErrorMap() {
        // POPULATE OUTPUT STACKS

        /// intensity-scaled stack
        FloatProcessor fpSRIntensityScaled = (FloatProcessor) fpSR.duplicate();
        fpSRIntensityScaled.multiply(alpha);
        fpSRIntensityScaled.add(beta);

        /// intensity-scaled and convolved stack
        FloatProcessor fpSRIntensityScaledBlurred = (FloatProcessor) fpSRIntensityScaled.duplicate();
        fpSRIntensityScaledBlurred.blurGaussian(sigma_linear);

        // CALCULATE METRICS AND MAP
        /// metrics
        float[] pixelsIntensityScaledBlurred_RefSize = (float[]) fpSRIntensityScaledBlurred.getPixels();
        this.globalRMSE = sqrt(calculateMSE(pixelsIntensityScaledBlurred_RefSize, pixelsRef));
        this.globalPPMCC = calculatePPMCC(pixelsIntensityScaledBlurred_RefSize, pixelsRef, true);

        /// error map
        float[] pixelsEMap = new float[width * height];
        float[] pixelsSRC = (float[]) fpSRIntensityScaledBlurred.getPixels();

        float maxRef = -Float.MAX_VALUE;
        for (int p = 0; p < width * height; p++) {
            float vRef = pixelsRef[p];
            float vSRC = pixelsSRC[p];

            maxRef = max(maxRef, vRef); //why
            pixelsEMap[p] = abs(vRef - vSRC);
        }

        return new FloatProcessor(width, height, pixelsEMap);


    }


    //////////////////////////////////////////////////////////////
    ///////////////// Helper functions & classes /////////////////
    //////////////////////////////////////////////////////////////


    private float[] calculateAlphaBeta(float[] xA, float[] y, float[] oneA) {
         /*
            xA = scaled and translated super-resolution
            oneA = scaled ones
            y =  reference
            */

        float N = 0;
        for (int i = 0; i < oneA.length; i++) {
            N += oneA[i] * oneA[i];
        }

        float nPixels = xA.length;
        assert (nPixels == y.length);
        assert (nPixels == oneA.length);


        float xATxA = 0, xAT1A = 0, yTxA = 0, yT1A = 0;

        for (int i = 0; i < nPixels; i++) {
            yTxA += y[i] * xA[i];
            yT1A += y[i] * oneA[i];
            xAT1A += xA[i] * oneA[i];
            xATxA += xA[i] * xA[i];
        }

        float numerator = N * yTxA - yT1A * xAT1A;
        float denominator = N * xATxA - xAT1A * xAT1A;
        float alphaHat = numerator / denominator;
        float betaHat = yT1A / N - alphaHat * (xAT1A / N);

        return new float[]{alphaHat, betaHat};
    }

    private float calculateRMSE(float[] array1, float[] array2) {

        int N = array1.length;
        double MSE = 0;

        for (int i = 0; i < N; i++) {
            MSE += (array1[i] - array2[i]) * (array1[i] - array2[i]);
        }
        MSE /= N;

        return (float) Math.sqrt(MSE);
    }


    private class sigmaOptimiseFunction implements UnivariateFunction {

        FloatProcessor fpSR;
        float[] pixelsRef, ones;
        ArrayList<Float> sigmaList = new ArrayList<Float>();
        ArrayList<Float> errorList = new ArrayList<Float>();

        public sigmaOptimiseFunction(FloatProcessor fpSR, float[] pixelsRef, float[] ones) {
            this.fpSR = fpSR;
            this.pixelsRef = pixelsRef;
            this.ones = ones;
        }

        public double value(double sigma) {
            FloatProcessor blurredFp = (FloatProcessor) fpSR.duplicate();
            FloatProcessor blurredOnes = new FloatProcessor(width, height, ones);
            blurredFp.blurGaussian(sigma);
            blurredOnes.blurGaussian(sigma);

//            blurredFp = (FloatProcessor) blurredFp.resize(w_Ref, h_Ref);
//            blurredOnes = (FloatProcessor) blurredOnes.resize(w_Ref, h_Ref);

            float[] aB = calculateAlphaBeta((float[]) blurredFp.getPixels(), pixelsRef, (float[]) blurredOnes.getPixels());

            FloatProcessor finalFpSR = (FloatProcessor) fpSR.duplicate();
            finalFpSR.multiply(aB[0]);
            finalFpSR.add(aB[1]);
            finalFpSR.blurGaussian(sigma);

//            FloatProcessor finalFpSRResized = (FloatProcessor) finalFpSR.resize(w_Ref, h_Ref);
//            double error = calculateRMSE(pixelsRef, (float[]) finalFpSRResized.getPixels());

            double error = calculateRMSE(pixelsRef, (float[]) finalFpSR.getPixels());

//            IJ.log("Optimising... sigma=" + df.format(sigma) + ", alpha=" + df.format(aB[0]) + ", beta=" + df.format(aB[1]) + ". Error=" + df.format(error));
            sigmaList.add((float) sigma);
            errorList.add((float) error);
            return error;
        }

        public ArrayList<Float> getSigmaList() {
            return sigmaList;
        }

        public ArrayList<Float> getErrorList() {
            return errorList;
        }
    }


    float getIntegratedGaussian(float dx, float dy, float sigma2) {
        float Ex = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
        float Ey = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
        float vKernel = Ex * Ey;
        return vKernel;
    }

    float erf(float g) {
        float x = abs(g);
        if (x >= 4.0f)
            return (g > 0.0f) ? 1.0f : -1.0f;

        // constants
        float a1 = 0.254829592f;
        float a2 = -0.284496736f;
        float a3 = 1.421413741f;
        float a4 = -1.453152027f;
        float a5 = 1.061405429f;
        float p = 0.3275911f;

        // A&S formula 7.1.26
        float t = 1.0f / (1.0f + p * x);
        float y = (float) (1.0f - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x));

        return (g > 0.0f) ? y : -y;
    }

    private FloatProcessor getRSF(float sigma) {

        // calculate the final RSF
        float sigma2 = (float) (ROOT2 * abs(sigma));
        int radius = max(((int) sigma) * 3, 1);
        int size = radius * 2 + 1;
        float vKernelSum = 0;

        FloatProcessor fpRSF = new FloatProcessor(size, size);
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                float vKernel = getIntegratedGaussian(dx, dy, sigma2);
                vKernelSum += vKernel;
                fpRSF.setf(dx + radius, dy + radius, vKernel);
            }
        }
        fpRSF.multiply(1. / vKernelSum);
        return fpRSF;
    }


}
