package nanoj.liveSRRF;

import ij.plugin.filter.GaussianBlur;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;


import static nanoj.core.java.image.calculator.FloatProcessorCalculator.multiply;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.subtract;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.add;
import static nanoj.core.java.image.calculator.FloatProcessorCalculator.divide;


public class SSIMCalculator {

    private FloatProcessor ipRef, ipRefGauss, ipRefSig2, ipRefGauss2;
    private GaussianBlur gaussBlurringMachine;
    private double radius;
    private final double accuracy = 0.01;
    private double C1, C2;


    // Initializor
    public SSIMCalculator(ImageProcessor ipRef, double radius, float[] regFactors){

        this.ipRef = ipRef.duplicate().convertToFloatProcessor();
        this.radius = radius;
        this.gaussBlurringMachine = new GaussianBlur();
        this.C1 = regFactors[0];
        this.C2 = regFactors[1];

        //Get mu_ref
        ipRefGauss = ipRef.duplicate().convertToFloatProcessor();
        gaussBlurringMachine.blurFloat(ipRefGauss, radius, radius, accuracy);

        // Get mu_ref^2
        ipRefGauss2 = multiply(ipRefGauss, ipRefGauss);

        // Get sigma_ref2
        FloatProcessor ipRef2Gauss = multiply(this.ipRef, this.ipRef);
        gaussBlurringMachine.blurFloat(ipRef2Gauss, radius, radius, accuracy);
        ipRefSig2 = subtract(ipRef2Gauss, multiply(ipRefGauss, ipRefGauss));

    }

    // Calculator
    public FloatProcessor Calculate(ImageProcessor ip){

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

}

