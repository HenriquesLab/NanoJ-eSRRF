package nanoj.liveSRRF;

import ij.ImageStack;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import static java.lang.Math.*;
import static nanoj.core.java.array.ArrayMath.calculateMSE;
import static nanoj.core.java.array.ArrayMath.calculatePPMCC;
import static nanoj.kernels.Kernel_BasePSO.*;
import static nanoj.kernels.Kernel_BasePSO.CONSTANT;

public class Build_RSE_RSP_maps {

    int maxMag = 4;
    int w_SR, h_SR, magnification;
    int currentFrame;

    double sigmaGuess, alphaGuess, betaGuess;
    float[] pixelsRef_scaledToSR;

    private boolean visualiseParameterEvolution = false;

    FloatProcessor fpSR, fpRef;
    private static Kernel_SquirrelSwarmOptimizer kPSOErrorMap = new Kernel_SquirrelSwarmOptimizer();
    ResultsTable rt = new ResultsTable();


    // --- Constructor ---
    public Build_RSE_RSP_maps(FloatProcessor fpSR, FloatProcessor fpRef, int w_SR, int h_SR, int magnification) {
        this.fpSR = fpSR;
        this.fpRef = fpRef;
        this.h_SR = h_SR;
        this.w_SR = w_SR;
        this.magnification = magnification;

        FloatProcessor fpRef_scaledToSR = (FloatProcessor) fpRef.duplicate();
        fpRef_scaledToSR.setInterpolationMethod(ImageProcessor.BICUBIC);
        fpRef_scaledToSR = (FloatProcessor) fpRef_scaledToSR.resize(w_SR, h_SR);
        this.pixelsRef_scaledToSR = (float[]) fpRef_scaledToSR.getPixels();

        this.currentFrame = 0;

    }

    // --- Calculate ---
    public ImageProcessor calculate() {

        int nPixelsSR = w_SR * h_SR;
//        ImageStack imsEMap = new ImageStack(w_SR, h_SR,1);

        // Set up SR float processors
        fpSR.resetRoi();
        FloatProcessor fpSR_scaledToRef = (FloatProcessor) fpSR.duplicate();
        fpSR_scaledToRef = (FloatProcessor) fpSR_scaledToRef.resize(fpRef.getWidth(), fpRef.getHeight());


//        log.msg("Meta-optimisation of alpha and beta:");
//
//        // Rough estimate of alpha and beta on reference image magnification
//        log.status("Meta-optimising alpha and beta...");


        CurveFitter cf = new CurveFitter(floatArrayToDoubleArray(fpSR_scaledToRef.getPixels()), floatArrayToDoubleArray(fpRef.getPixels()));
        cf.doFit(cf.STRAIGHT_LINE);
        betaGuess = cf.getParams()[0];
        alphaGuess = cf.getParams()[1];
//        log.msg("\t Initial alpha guess = " + alphaGuess + ", initial beta guess = " + betaGuess);

//        if (impRSF != null)
//
//        {
//            log.status("Extracting RSF features...");
//            FloatProcessor fpRSF = imsRSF.getProcessor(min(n, nSlicesRSF)).convertToFloatProcessor();
//        GaussianFitMinimizer gaussianFitMinimizer = new GaussianFitMinimizer(fpRSF, 1.75, fpRSF.getWidth() / 2, fpRSF.getHeight() / 2);
//        Double[] fitResults = gaussianFitMinimizer.calculate();
//        sigmaGuess = fitResults[2];
//        IJ.log("Extracted sigma from fit to RSF is: " + sigmaGuess);
//        } else

//        {
        sigmaGuess = 0;
//        }

        // Main optimizer - estimate alpha, beta and sigma

//        log.msg("Joint PSO Optimisation:");
//        log.status("Joint optimisation of alpha, beta and sigma...");
        double[] results;

        kPSOErrorMap.maxMagnification = maxMag;

        if (alphaGuess != 0 && betaGuess != 0)

        {
            results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
                    new double[]{0.001, betaGuess * 0.1, 0.5}, // low boundary
                    new double[]{alphaGuess * 10, betaGuess * 10, 10}, // high boundary
                    new double[]{alphaGuess, betaGuess, sigmaGuess == 0 ? 5 : sigmaGuess}, // best guess
                    new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, sigmaGuess == 0 ? OBEY_BOUNDARY : CONSTANT}, // impose boundaries
                    1e-7);
        } else

        {
            results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
                    new double[]{0.001, -1000, 0.5}, // low boundary
                    new double[]{100, 1000, 10}, // high boundary
                    new double[]{10, 0, sigmaGuess == 0 ? 5 : sigmaGuess}, // best guess
                    new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, sigmaGuess == 0 ? OBEY_BOUNDARY : CONSTANT}, // impose boundaries
                    Double.NaN);
        }

//            if (sigmaGuess == 0) {
//                if (alphaGuess != 0 && betaGuess != 0) {
//                    results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
//                            new double[]{0.001, betaGuess * 0.1, 0.5}, // low boundary
//                            new double[]{alphaGuess * 10, betaGuess * 10, 10}, // high boundary
//                            new double[]{alphaGuess, betaGuess, 5}, // best guess
//                            new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, OBEY_BOUNDARY}, // impose boundaries
//                            1e-7);
//                }
//                else {
//                    results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
//                            new double[]{0.001, -1000, 0.5}, // low boundary
//                            new double[]{100, 1000, 10}, // high boundary
//                            new double[]{10, 0, 5}, // best guess
//                            new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, OBEY_BOUNDARY}, // impose boundaries
//                            Double.NaN);
//                }
//            }
//            else {
//                if (alphaGuess != 0 && betaGuess != 0) {
//                    results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
//                            new double[]{0.001, betaGuess * 0.1, sigmaGuess * 0.9}, // low boundary
//                            new double[]{alphaGuess * 10, betaGuess * 10, sigmaGuess * 1.1}, // high boundary
//                            new double[]{alphaGuess, betaGuess, sigmaGuess}, // best guess
//                            new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, OBEY_BOUNDARY}, // impose boundaries
//                            1e-7);
//                } else {
//                    results = kPSOErrorMap.calculate(fpRef, fpSR, 10, 300,
//                            new double[]{0.001, -1000, sigmaGuess * 0.9}, // low boundary
//                            new double[]{100, 1000, sigmaGuess * 1.1}, // high boundary
//                            new double[]{10, 0, sigmaGuess}, // best guess
//                            new int[]{OBEY_LOW_BOUNDARY, DONT_OBEY_BOUNDARY, OBEY_BOUNDARY}, // impose boundaries
//                            Double.NaN);
//                }
//            }

        double alpha = results[0];
        double beta = results[1];
        double sigma = results[2];

        if (visualiseParameterEvolution)

        {
            // Visualization of parameter evolution
            kPSOErrorMap.plotOptimizationEvolution(new String[]{"Alpha", "Beta", "Sigma"});
            kPSOErrorMap.plotOptimizationEvolution3D(512, 512, 0, 1, 2).show();
            kPSOErrorMap.renderOptimizationEvolution2D(128, 128, 0, 1).show();
            kPSOErrorMap.renderOptimizationEvolution2D(128, 128, 0, 2).show();
        }

//        log.msg("Final parameters:");
//        log.msg("\t Alpha = " + alpha + ", Beta = " + beta + ", Sigma = " + sigma + " (in Ref pixels)");

        // Generate RSF for this image
//        log.status("Generating RSF");
        FloatProcessor fpRSF = kPSOErrorMap.getRSF();
//        fpRSFArray[n - 1] = fpRSF;
//        maxRSFWidth =
//
//                max(maxRSFWidth, fpRSF.getWidth());

//        log.status("Intensity rescaling SR image");
        float[] pixelsSR_intensityMatched = (float[]) fpSR.getPixels();
        for (
                int p = 0;
                p < nPixelsSR; p++)

        {
            pixelsSR_intensityMatched[p] = (float) (pixelsSR_intensityMatched[p] * alpha + beta);
        }

//        imsSRNormalised.setProcessor(new FloatProcessor(w_SR, h_SR, pixelsSR_intensityMatched), n);

        // Populate intensity rescaled and convolved SR image stack
//        log.status("Convolving RSF with SR");
        FloatProcessor fpSRC = (FloatProcessor) fpSR.duplicate();
        fpSRC.blurGaussian(sigma * magnification);

//        imsSRConvolved.setProcessor(fpSRC, n);

        // Put into error map

//        log.status("Calculating similarity...");

        int w_Ref = fpRef.getWidth();
        int h_Ref = fpRef.getHeight();

        FloatProcessor fpSRC_scaledToRef = (FloatProcessor) fpSRC.duplicate();
        fpSRC_scaledToRef = (FloatProcessor) fpSRC_scaledToRef.resize(w_Ref, h_Ref);
        float[] pixelsSRC_scaledToRef = (float[]) fpSRC_scaledToRef.getPixels();
        float[] pixelsRef = (float[]) fpRef.getPixels();

        double globalRMSE = sqrt(calculateMSE(pixelsSRC_scaledToRef, pixelsRef));
        double globalPPMCC = calculatePPMCC(pixelsSRC_scaledToRef, pixelsRef, true);

        float[] pixelsEMap = new float[nPixelsSR];
        float[] pixelsSRC = (float[]) fpSRC.getPixels();

        float maxRef = -Float.MAX_VALUE;
        for (
                int p = 0;
                p < pixelsEMap.length; p++)

        {
            float vRef = pixelsRef_scaledToSR[p];
            float vSRC = pixelsSRC[p];

            maxRef = max(maxRef, vRef);
//            if (showPositiveNegative) pixelsEMap[p] = (vRef - vSRC);
//            else pixelsEMap[p] = abs(vRef - vSRC);
            pixelsEMap[p] = abs(vRef - vSRC);
        }

        // set Error Map into stack
//        imsEMap.setProcessor(new FloatProcessor(w_SR, h_SR, pixelsEMap), 1);

        FloatProcessor fpErrorMap = new FloatProcessor(w_SR, h_SR, pixelsEMap);

        // present error in table
        rt.incrementCounter();
        rt.addValue("Frame", currentFrame);
        rt.addValue("RSP (Resolution Scaled Pearson-Correlation)", globalPPMCC);
        rt.addValue("RSE (Resolution Scaled Error)", globalRMSE);

//        if (nSlicesSR > 1)
//
//        {
//            double frameTime = ((System.nanoTime() - loopStart) / n) / 1e9;
//            double remainingTime = frameTime * (nSlicesSR - n);
//            int _h = (int) (remainingTime / 3600);
//            int _m = (int) (((remainingTime % 86400) % 3600) / 60);
//            int _s = (int) (((remainingTime % 86400) % 3600) % 60);
////            log.msg("Estimated time remaining to complete analysis:");
////            log.msg("\t" + String.format("%02d:%02d:%02d", _h, _m, _s));
//        }

        currentFrame++;

//        return imsEMap.getProcessor(1);
        return fpErrorMap;

    }


    private double[] floatArrayToDoubleArray(Object pixels) {
        float[] floatArray = (float[]) pixels;
        double[] doubleArray = new double[floatArray.length];
        for (int i = 0; i < floatArray.length; i++) {
            doubleArray[i] = (double) floatArray[i];
        }
        return doubleArray;
    }


}
