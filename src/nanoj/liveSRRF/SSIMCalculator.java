package nanoj.liveSRRF;

import ij.IJ;
import ij.ImageStack;
import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


import java.util.Arrays;

import static nanoj.core.java.image.calculator.FloatProcessorCalculator.multiply;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.subtract;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.divide;
import static nanoj.core2.NanoJImageStackArrayConvertion.ImageStackToFloatArray;


public class SSIMCalculator {

    private FloatProcessor ipRef, ipRefGauss, ipRefSig2, ipRefGauss2;
    private GaussianBlur gaussBlurringMachine;
    private double radius;
    private final double accuracy = 0.01;
    private final double C1, C2;
//    private double mu_ref, sig_ref;
    private double[][] musig;
    private ImageStack ims;

    private static boolean DEBUG = false;


    // Constructor #1
    public SSIMCalculator(ImageProcessor ipRef, double radius, float[] regFactors){

        this.ipRef = ipRef.duplicate().convertToFloatProcessor();
        this.radius = radius;
        this.C1 = regFactors[0];
        this.C2 = regFactors[1];

        if (radius < Float.MAX_VALUE) {
            //Get mu_ref
            if (DEBUG) {
                IJ.log("Initialising SSIM calculator...");
            }
            this.gaussBlurringMachine = new GaussianBlur();
            ipRefGauss = ipRef.duplicate().convertToFloatProcessor();
            gaussBlurringMachine.blurFloat(ipRefGauss, radius, radius, accuracy);

            // Get mu_ref^2
            ipRefGauss2 = multiply(ipRefGauss, ipRefGauss);

            // Get sigma_ref2
            FloatProcessor ipRef2Gauss = multiply(this.ipRef, this.ipRef);
            gaussBlurringMachine.blurFloat(ipRef2Gauss, radius, radius, accuracy);
            ipRefSig2 = subtract(ipRef2Gauss, multiply(ipRefGauss, ipRefGauss));
        }

//        double[] musig_ref = getMeanStdofImage(ipRef);
//        mu_ref = musig_ref[0];
//        sig_ref = musig_ref[1];
    }


    // Constructor #2
    public SSIMCalculator(ImageStack ims, float[] regFactors){

        this.C1 = regFactors[0];
        this.C2 = regFactors[1];
        this.ims = ims;

        musig = new double[ims.size()][2];

        for (int i = 0; i < ims.size(); i++) {
            musig[i] = getMeanStdofImage(ims.getProcessor(i+1).convertToFloatProcessor());
        }
    }


//    public double CalculateMetric(ImageProcessor ip){
//        double[] musig = getMeanStdofImage(ip);
//        double[] musigCC = getMeanStdofImage(multiply(ip.duplicate().convertToFloatProcessor(), ipRef));
//        double CC = musigCC[0]-musig[0]*mu_ref;
//
//        double ssim = ((2.0f*mu_ref*musig[0] + C1)*(2.0f*CC + C2)/((mu_ref*mu_ref + musig[0]*musig[0] + C1)*(sig_ref*sig_ref + musig[1]*musig[1] + C2)));
//
//        return ssim;
//    }

    public double[] rollingSSIM(int refindex, int traceLength){

        FloatProcessor fpRef = ims.getProcessor(refindex+1).duplicate().convertToFloatProcessor();
        double[] ssimTrace = new double[traceLength];
        double[] musigCC;
        double CC;
        double muref2 = musig[refindex][0]*musig[refindex][0];
        double musig2 = musig[refindex][1]*musig[refindex][1];
        for (int i = 0; i < traceLength; i++) {
            musigCC = getMeanStdofImage(multiply(ims.getProcessor(refindex+i+2).duplicate().convertToFloatProcessor(), fpRef));
            CC = musigCC[0]-musig[refindex+1+i][0]*musig[refindex][0];
            ssimTrace[i] = ((2.0f*musig[refindex][0]*musig[refindex+1+i][0] + C1)*(2.0f*CC + C2)/((muref2 + musig[refindex+1+i][0]*musig[refindex+1+i][0] + C1)*(musig2 + musig[refindex+1+i][1]*musig[refindex+1+i][1] + C2)));
        }
        return ssimTrace;
    }


//    public double[] CalculateMetric(ImageStack ims){
//        // TODO: This doesn't work at the moment
//        NanoJThreadExecutor NTE = new NanoJThreadExecutor(false);
//
//        double[] ssimTrace = new double[ims.size()];
//        for (int i = 0; i < ims.size(); i++) {
////            ssimTrace[i] = CalculateMetric(ims.getProcessor(i+1).convertToFloatProcessor());
//            ThreadedSSIM t = new ThreadedSSIM(ims, ssimTrace, i+1);
//            NTE.execute(t);
//        }
//
//        NTE.finish();
//        return ssimTrace;
//    }



    // Calculator
    public FloatProcessor CalculateSSIMmap(ImageProcessor ip){

        // Get mu
        FloatProcessor ipGauss = ip.duplicate().convertToFloatProcessor();
        gaussBlurringMachine.blurFloat(ipGauss, radius, radius, accuracy);

        // Get mu^2
        FloatProcessor ipGauss2 = multiply(ipGauss, ipGauss);

        // Get sigma2
        FloatProcessor ip2Gauss = multiply(ip.duplicate().convertToFloatProcessor(), ip.duplicate().convertToFloatProcessor());
        gaussBlurringMachine.blurFloat(ip2Gauss, radius, radius, accuracy);
        FloatProcessor ipSig2 = subtract(ip2Gauss, multiply(ipGauss, ipGauss));

        // Get cross-covariance
        FloatProcessor ipipRefGauss = multiply(ip.duplicate().convertToFloatProcessor(), ipRef);
        gaussBlurringMachine.blurFloat(ipipRefGauss, radius, radius, accuracy);

        FloatProcessor cc = subtract(ipipRefGauss, multiply(ipRefGauss, ipGauss));

        // Calculate SSIM
        FloatProcessor A = multiply(ipRefGauss, ipGauss);
        A.multiply(2.0f);
        A.add(C1);

        FloatProcessor B = cc.duplicate().convertToFloatProcessor();
        B.multiply(2.0f);
        B.add(C2);

        FloatProcessor C = add(ipRefGauss2, ipGauss2);
        C.add(C1);

        FloatProcessor D = add(ipSig2, ipRefSig2);
        D.add(C2);

        FloatProcessor ssim = divide(multiply(A, B), multiply(C, D));

        return ssim;
    }

    // ------------------ Functions ------------------
    public static float[] getRegularizationFactors(ImageStack ims, String regMethod){

        float[] C = new float[2];
        // Regularization factors
        if (regMethod.equals("image bit depth")) {
            int bitDepth = ims.getBitDepth();
            C[0] = (float) Math.pow(0.01f * (float) (Math.pow(2.0f, bitDepth) - 1), 2);
            if (DEBUG) {
                IJ.log("Bit depth: " + bitDepth);
            }
        }
        else{
            float[] percentiles = new float[2];
            percentiles[0] = getPercentileFromImageStack(ims, 0.01f);
            percentiles[1] = getPercentileFromImageStack(ims, 0.99f);
            if (DEBUG) {
                IJ.log("1% percentile: " + percentiles[0]);
                IJ.log("99% percentile: " + percentiles[1]);
            }
            C[0] = (float) Math.pow(0.01f * (percentiles[1] - percentiles[0]), 2);
        }
        C[1] = 9 * C[0];

        if (DEBUG) {
            IJ.log("Regularization factors: " + regMethod);
            IJ.log("C1: " + C[0]);
            IJ.log("C2: " + C[1]);
        }

        return C;

    }

    public static float getPercentileFromImageStack(ImageStack ims, float percentile){

        float[] pixelValues = ImageStackToFloatArray(ims);
        Arrays.sort(pixelValues);
        int nPixels = ims.getWidth()*ims.getHeight()*ims.getSize();
        int id = (int) ((float) nPixels*percentile);

        return pixelValues[id];
    }

    public static double[] getMeanStdofImage(ImageProcessor ip){

        float[] values = (float[]) ip.duplicate().getPixels();
        double m = 0;
        double s = 0;
        for (float value : values) {
            m += value;
            s += value * value;
        }
        m /= values.length;
        s /= values.length;
        s = Math.sqrt(s - m*m);

        return new double[]{m, s};
    }
//    class ThreadedSSIM extends Thread {
//        private final ImageProcessor ip;
//        private final double[] ssim;
//        private final int n;
//
//
//        public ThreadedSSIM(ImageStack ims, double[] ssimArray, int n) {
//            this.n = n;
//            this.ssim = ssimArray;
//            this.ip = ims.getProcessor(n).convertToFloatProcessor();
//        }
//
//        @Override
//        public void run() {
//            float[] musig = getMeanStdofImage(ip);
//            float[] musigCC = getMeanStdofImage(multiply(ip.duplicate().convertToFloatProcessor(), ipRef.convertToFloatProcessor()));
//            float CC = musigCC[0]-musig[0]*mu_ref;
//
//            ssim[n] = ((2.0f*mu_ref*musig[0] + C1)*(2.0f*CC + C2)/((mu_ref*mu_ref + musig[0]*musig[0] + C1)*(sig_ref*sig_ref + musig[1]*musig[1] + C2)));
//
//        }
//    }



}

