package nanoj.liveSRRF;

import ij.measure.Minimizer;
import ij.measure.UserFunction;
import ij.process.FloatProcessor;

import static java.lang.Math.*;

/**
 * Created by paxcalpt on 01/06/2017.
 */
public class GaussianFitMinimizer implements UserFunction {

    private final FloatProcessor ip;
    private final double initialSigma;
    private final double initialX;
    private final double initialY;
    private final double initialBackground;
    private final double initialAmplitude;

    private final double maxSigma;
    private final int w, h;
    public double sigma, background, amplitude, xc, yc, error;

    public GaussianFitMinimizer(FloatProcessor ip, double initialSigma, double initialX, double initialY) {
        this.ip = ip;
        this.w = ip.getWidth();
        this.h = ip.getHeight();
        this.initialSigma = initialSigma;
        this.initialX = initialX;
        this.initialY = initialY;
        this.initialBackground = ip.getMin();
        this.initialAmplitude = ip.getMax() - initialBackground;
        this.maxSigma = (min(w,h)-1d)/2d;
    }

    public Double[] calculate() {
        double[] initialParameters = new double[] {initialSigma, initialX, initialY, initialAmplitude, initialBackground};
        double[] initialParametersVariation = new double[] {0.1, 0.5, 0.5, 10, 10};

        Minimizer min = new Minimizer();
        min.setFunction(this, 5);
        min.setMaxError(0);
        //min.setStatusAndEsc("Estimating PSF: Iteration ", true);
        min.minimize(initialParameters, initialParametersVariation);
        sigma = min.getParams()[0];
        xc = min.getParams()[1];
        yc = min.getParams()[2];
        amplitude = min.getParams()[3];
        background = min.getParams()[4];
        error = calculateError(xc, yc, 2*sigma*sigma, amplitude, background);

        return new Double[] {xc, yc, sigma, amplitude, background, error};
    }

    @Override
    public double userFunction(double[] params, double p) {

        //System.out.println(Arrays.toString(params));
        double sigma = params[0];
        double sigma22 = 2*sigma*sigma;
        double xc = params[1];
        double yc = params[2];
        double amplitude = params[3];
        double background = params[4];

        if(sigma<0) return Double.NaN;
        if(sigma>maxSigma) return Double.NaN;
        if(xc<0.25*w || yc<0.25*h) return Double.NaN;
        if(xc>0.75*w || yc>0.75*h) return Double.NaN;

        return calculateError(xc, yc, sigma22, amplitude, background);
    }

    private double calculateError(double xc, double yc, double sigma22, double amplitude, double background) {
        double error = 0;
        for (int j=0; j<h; j++) {
            for (int i=0; i<w; i++) {
                double gv = amplitude*exp(-(pow(i-xc,2)+pow(j-yc,2))/sigma22)+background;
                double v = ip.getf(i, j);
                error += pow(gv-v, 2);
            }
        }
        return error;
    }
}
