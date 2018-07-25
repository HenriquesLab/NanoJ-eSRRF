package nanoj.liveSRRF;

import ij.IJ;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import nanoj.kernels.Kernel_BasePSO;

import java.util.Random;

/**
 * Created by Henriques-lab on 25/06/2017.
 */
public class Kernel_SquirrelSwarmOptimizer extends Kernel_BasePSO {

    private float ROOT2 = sqrt(2);

    private static final int STEP_CONVOLVE_HORIZONTAL = 0;
    private static final int STEP_CONVOLVE_VERTICAL = 1;
    private static final int STEP_CALCULATE_ERROR_MAP = 2;
    private int step = STEP_CONVOLVE_HORIZONTAL;

    private int widthRef, heightRef, widthHeightRef, widthSR, heightSR, widthHeightSR, magnification, magnification2, originalMagnification;

    private float[] pixelsRef, pixelsSR;
    private float[] particlesErrorMap, particlesSRC_H, particlesSRC_V;

    public int maxMagnification = 5;

    // Note: swarmSize can't be too high!! this uses > (fpRef + fpSR + fpSR*swarmSize) memory
    public double[] calculate(FloatProcessor fpRef, FloatProcessor fpSR, int swarmSize, int nIterations,
                              double[] position_lowBoundary, double[] position_highBoundary, double[] position_bestGuess, int[] position_constrainToBoundaries, double minimalImprovement) {

        this.originalMagnification = fpSR.getWidth() / fpRef.getWidth();
        if (this.originalMagnification > maxMagnification) { // too high magnification, truncating it to 5
            this.magnification = maxMagnification;
            fpSR.setInterpolationMethod(ImageProcessor.BICUBIC);
            fpSR = (FloatProcessor) fpSR.resize(fpRef.getWidth()*maxMagnification, fpRef.getHeight()*maxMagnification);
        }
        else
            this.magnification = this.originalMagnification;

        this.generator = new Random(327636000); // always use same seed for repeatability;
        this.widthRef = fpRef.getWidth();
        this.heightRef = fpRef.getHeight();
        this.widthHeightRef = widthRef * heightRef;
        double aspectRatioRef = widthRef / heightRef;
        this.widthSR = fpSR.getWidth();
        this.heightSR = fpSR.getHeight();
        this.widthHeightSR = widthSR * heightSR;
        double aspectRatioSR = widthSR / heightSR;
        this.magnification = widthSR / widthRef;
        this.magnification2 = this.magnification * this.magnification;

        assert (aspectRatioRef == aspectRatioSR);
        assert (position_lowBoundary.length == position_highBoundary.length);
        assert (position_bestGuess.length == position_lowBoundary.length);

        this.widthRef = fpRef.getWidth();
        this.heightRef = fpRef.getHeight();
        this.widthHeightRef = widthRef * heightRef;

        // initialise the PSO engine
        this.initializePSO(swarmSize, nIterations,
                position_lowBoundary, position_highBoundary,
                position_bestGuess, position_constrainToBoundaries,
                minimalImprovement);

        // Initialise Arrays
        this.pixelsRef = (float[]) fpRef.getPixels();
        this.pixelsSR = (float[]) fpSR.getPixels();
        this.particlesSRC_H = new float[swarmSize * widthHeightSR];
        this.particlesSRC_V = new float[swarmSize * widthHeightSR];
        this.particlesErrorMap = new float[swarmSize * widthHeightRef];

        // Upload arrays
        setExplicit(true);
        autoChooseDeviceForNanoJ();
        put(this.pixelsRef);
        put(this.pixelsSR);
        put(this.particlesSRC_H);
        put(this.particlesSRC_V);
        put(this.particlesErrorMap);

        // Implement the PSO

        long timeStart = System.nanoTime();
        long timeLastUpdate = System.currentTimeMillis();

        for (int iteration = 1; iteration <= nIterations; iteration++) {
            if (this.stop || IJ.escapePressed()) {
                IJ.resetEscape();
                log.status("done...");
                log.progress(1);
                this.nIterations = iteration-1;
                break;
            }

            log.progress(iteration, nIterations);

            // Calculate Error-Map on GPU
            put(this.position_particles); // note that in each new cycle we re-upload the new calculated positions
            step = STEP_CONVOLVE_HORIZONTAL;
            executeByBlocks(swarmSize * widthHeightSR); // !! EXECUTE OPENCL CODE !!
            //execute(swarmSize * widthHeightSR); // !! EXECUTE OPENCL CODE !!
            step = STEP_CONVOLVE_VERTICAL;
            executeByBlocks(swarmSize * widthHeightSR); // !! EXECUTE OPENCL CODE !!
            //execute(swarmSize * widthHeightSR); // !! EXECUTE OPENCL CODE !!
            step = STEP_CALCULATE_ERROR_MAP;
            executeByBlocks(swarmSize * widthHeightRef); // !! EXECUTE OPENCL CODE !!
            //execute(swarmSize * widthHeightRef); // !! EXECUTE OPENCL CODE !!
            get(this.particlesErrorMap); // pull the error map from the GPU

            // Calculate Error and Update Particles Positions

            for(int pId=0; pId<swarmSize; pId++) {

                // Calculate the error
                double error = 0;
                int counter = 1;
                for (int n = pId* widthHeightRef; n<(pId+1)* widthHeightRef; n++) {
                    double vError = this.particlesErrorMap[n];
                    error += (vError-error) / counter;
                    counter++;
                }
                error = sqrt(error);

                // Update Global and Local Error Values
                if (!Double.isNaN(error)) { // make sure error is not funky
                    this.updateError(error, iteration, pId);
                }
            }

            // Update Particles Position
            this.updateParticlePosition(iteration);

            if (System.currentTimeMillis() - timeLastUpdate > 1000) { // update every second
                timeLastUpdate = System.currentTimeMillis();
                double iterationTime= ((System.nanoTime()-timeStart)/(iteration+1))/1e9;
                double remainingTime = iterationTime*(nIterations-iteration);
                int _h = (int) (remainingTime / 3600);
                int _m = (int) (((remainingTime % 86400) % 3600) / 60);
                int _s = (int) (((remainingTime % 86400) % 3600) % 60);
                String abs = String.format("a=%.2f b=%.2f s=%.2f ",
                        globalBestPosition[0], globalBestPosition[1], globalBestPosition[2]);
                String etf;
                if (_h > 0) etf = String.format("ETF %02d:%02d:%02d ", _h, _m, _s);
                else if (_m > 0) etf = String.format("ETF %02d:%02d ", _m, _s);
                else etf = String.format("ETF %ds ", _s);
                log.status(abs+etf+String.format(" %.1f", 100. * iteration / nIterations)+"%");
            }
        }

        return globalBestPosition;
    }

    @Override
    public void run() {
        if (step == STEP_CONVOLVE_HORIZONTAL) convolveHorizontal();
        else if (step == STEP_CONVOLVE_VERTICAL) convolveVertical();
        else if (step == STEP_CALCULATE_ERROR_MAP) calculateErrorMap();
    }

    public void convolveHorizontal() {
        int p = getGlobalId(0)+this.blockOffset;
        int xSR = p % widthSR;
        int ySR = (p / widthSR) % heightSR;
        int pId = p / (widthHeightSR);

        int p0v = pId * nVariables;
        float alpha = position_particles[p0v];
        float beta  = position_particles[p0v + 1];
        float sigma = position_particles[p0v + 2] * magnification;
        float sigma2 = ROOT2*abs(sigma); // 2 * pow(sigma, 2);

        int radius = max(((int) sigma) * 3, 1);

        float vKernelSum = 0;
        float v = 0;

        for (int dx = -radius; dx <= radius; dx++) {
            float vKernel = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
            int xp = min(max(xSR+dx, 0), widthSR-1);
            v += (pixelsSR[ySR * widthSR + xp] * alpha + beta) * vKernel;
            vKernelSum += vKernel;
        }

        v /= vKernelSum;
        this.particlesSRC_H[p] = v;
    }

    public void convolveVertical() {
        int p = getGlobalId(0)+this.blockOffset;
        int xSR = p % widthSR;
        int ySR = (p / widthSR) % heightSR;
        int pId = p / (widthHeightSR);

        int p0v = pId * nVariables;
        float sigma = position_particles[p0v + 2] * magnification;
        float sigma2 = ROOT2*abs(sigma); // 2 * pow(sigma, 2);

        int radius = max(((int) sigma) * 3, 1);
        int offset = widthHeightSR * pId;

        float vKernelSum = 0;
        float v = 0;

        for (int dy = -radius; dy <= radius; dy++) {
            float vKernel = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
            int yp = min(max(ySR+dy, 0), heightSR-1);
            v += particlesSRC_H[offset + yp * widthSR + xSR] * vKernel;
            vKernelSum += vKernel;
        }

        v /= vKernelSum;
        this.particlesSRC_V[p] = v;
    }

    public void calculateErrorMap() {
        int p = getGlobalId(0)+this.blockOffset;

        int xRef = p % widthRef;
        int yRef = (p / widthRef) % heightRef;
        int pId = p / (widthHeightRef);

        int offset = widthHeightSR * pId;

        float v = 0;
        for (int j = 0; j < magnification; j++) {
            for (int i = 0; i < magnification; i++) {
                int xSR = xRef*magnification + i;
                int ySR = yRef*magnification + j;
                v += particlesSRC_V[offset + widthSR * ySR + xSR];
            }
        }
        v /=  magnification2;
        particlesErrorMap[p] = pow(v - pixelsRef[yRef* widthRef +xRef],2);
    }

    public FloatProcessor getRSF() {
        // calculate the final RSF
        double sigma = globalBestPosition[2] * originalMagnification;
        float sigma2 = (float) (ROOT2*abs(sigma));
        int radius = max(((int) sigma) * 3, 1);
        int size = radius * 2 + 1;
        float vKernelSum = 0;

        FloatProcessor fpRSF = new FloatProcessor(size, size);
        for (int dy = -radius; dy <= radius; dy++) {
            for (int dx = -radius; dx <= radius; dx++) {
                float vKernel = getIntegratedGaussian(dx, dy, sigma2);
                vKernelSum+= vKernel;
                fpRSF.setf(dx+radius, dy+radius, vKernel);
            }
        }
        fpRSF.multiply(1./vKernelSum);
        return fpRSF;
    }

    // based on Smith et al, Nath Meth, 2010: Fast, single-molecule localization that achieves theoretically
    // minimum uncertainty (see Sup Mat page 10)
    // note, the paper has an error on their formula 4a and 4b, 2sigma^2 should be sqrt(2)*sigma
    // see https://en.wikipedia.org/wiki/Normal_distribution formula 'Cumulative distribution function'
    float getIntegratedGaussian(float dx, float dy, float sigma2) {
        float Ex = 0.5f * (erf((dx + 0.5f) / sigma2) - erf((dx - 0.5f) / sigma2));
        float Ey = 0.5f * (erf((dy + 0.5f) / sigma2) - erf((dy - 0.5f) / sigma2));
        float vKernel = Ex * Ey;
        return vKernel;
    }

    /*
    * The following code is adapted from
    * http://www.johndcook.com/blog/2009/01/19/stand-alone-error-function-erf/
    * which is licensed Public Domain.
    *
    * Limited precision, only useful for floats. Maximum error is below 1.5 × 10-7.
    */
    float erf(float g) {
        float x = abs(g);
        if (x >= 4.0f)
            return (g > 0.0f) ? 1.0f : -1.0f;

        // constants
        float a1 =  0.254829592f;
        float a2 = -0.284496736f;
        float a3 =  1.421413741f;
        float a4 = -1.453152027f;
        float a5 =  1.061405429f;
        float p  =  0.3275911f;

        // A&S formula 7.1.26
        float t = 1.0f / (1.0f + p*x);
        float y = 1.0f - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

        return (g > 0.0f) ? y : -y;
    }

}
