package nanoj.liveSRRF;

import ij.measure.Minimizer;
import ij.measure.UserFunction;
import org.python.modules.math;

public class Fit1DGaussian implements UserFunction {

    public float initAmp, initSigma, initBG, fittedAmp, fittedSigma, fittedBG, fittedX0;
    public int initX0, nElements, startCrop, endCrop;
    public float[] dataArray;

    // Initilization of the fitting
    public Fit1DGaussian(float[] array) {

        this.nElements = array.length;
        this.dataArray = new float[nElements];
        this.startCrop = 0;
        this.endCrop = nElements-1;


        this.initAmp = array[0];
        this.initX0 = 0;
        this.initBG = array[0];

        int xMin = 0;

        for (int x = 0; x < nElements; x++) {
            if (array[x] > initAmp) {
                initAmp = array[x];
                this.initX0 = x;
            }
            if (array[x] < initBG) {
                initBG = array[x];
                xMin = x;
            }
        }
        initAmp = initAmp - initBG;

        int initSigmaX = 0;
        if (xMin < initX0) while ((array[xMin + initSigmaX]-initBG)/initAmp < 0.6065) initSigmaX ++; // rising edge
        else while ((array[initX0 + initSigmaX]-initBG)/initAmp > 0.6065) initSigmaX ++; // descending edge exp(-0.5) ~ 0.6065
        this.initSigma = (float) math.fabs(initSigmaX-initX0);

        for (int x = 0; x < nElements; x++){
            dataArray[x] = (array[x] - initBG)/initAmp; // normalizing the array
        }

        this.fittedAmp = initAmp;
        this.fittedX0 = initX0;
        this.fittedSigma = initSigma;
        this.fittedBG = initBG;
    }

    // Crop the array to fit only the top of the curve
    public void cropDataArray(float CropLevel){
        startCrop = initX0;
        endCrop = initX0;
        while (dataArray[startCrop] > CropLevel && startCrop > 0) startCrop--;
        while (dataArray[endCrop] > CropLevel && endCrop < nElements-1) endCrop++;

    }

    // Carry out the fitting
    public float[] calculate() {

        double[] initialParameters = new double[] {1, initX0, initSigma, 0};
        double[] initialParametersVariation = new double[] {0.1, 1, 1, 0.1}; // something smarter might be possible

        Minimizer min = new Minimizer();
        min.setFunction(this, 4);
        min.setMaxIterations(1000);
        min.setStatusAndEsc("Estimating max subpixel position: Iteration ", true);
        min.minimize(initialParameters, initialParametersVariation);
        this.fittedAmp = (float) min.getParams()[0]*initAmp;
        this.fittedX0 = (float) min.getParams()[1];
        this.fittedSigma = (float) min.getParams()[2];
        this.fittedBG = (float) min.getParams()[3]*initAmp + initBG;

        return new float[]{fittedAmp, fittedX0, fittedSigma, fittedBG};
    }

    // Return a float array for plotting (double[] because it feeds directly into plot)
    public double[][] fittedCurve() {

        float sigma22 = (float) (2*math.pow(fittedSigma,2));
        double[][] fittedCurveArray = new double[2][nElements];

        for (int x = 0; x < nElements; x++) {
            fittedCurveArray[0][x] = (double) initAmp* dataArray[x] + initBG;
            fittedCurveArray[1][x] = fittedAmp*math.exp(-math.pow(x-fittedX0,2)/sigma22) + fittedBG;
        }
        return fittedCurveArray;
    }

    // This userFunction calculates the square residuals ---------------------------------------------------------------
    @Override
    public double userFunction(double[] params, double v) {
        double amp = params[0];
        double x0 = params[1];
        double sigma = params[2];
        double background = params[3];

        if (amp < 0) return Double.MAX_VALUE; // add penalty if amp goes beyond bounds, same could be applied to other parameters

        double sigma22 = 2*math.pow(sigma,2);

        float residuals = 0;
        for (int x = startCrop; x < endCrop +1; x++) {
            residuals += math.pow(dataArray[x] - amp* math.exp(-math.pow(x-x0, 2)/sigma22) - background, 2);
        }
        return residuals;
    }
}
