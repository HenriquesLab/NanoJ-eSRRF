package nanoj.liveSRRF;

import ij.IJ;
import ij.ImageStack;
import ij.gui.Plot;

import java.awt.*;

import static nanoj.core.java.array.ArrayMath.*;

public class GetAxialPositionMFM {

    private ImageStack imsRef;
    int blockSize, width, height;
    Plot defocusPlot;
    double[] zPosArray;
    Fit1DGaussian fitting;
    float zStep;

    public GetAxialPositionMFM(ImageStack imsRef, float zStep, int nSlices) {
        this.imsRef = imsRef;
        this.blockSize = imsRef.getSize();
        IJ.log("Size of reference block: "+imsRef.getSize());
        IJ.log("Size of total block: "+nSlices);
        this.width = imsRef.getWidth();
        this.height = imsRef.getHeight();
        this.defocusPlot = new Plot("Defocus plot", "Z (nm)", "Intensity (AU)");
        this.zPosArray = new double[nSlices-blockSize];
        for (int n = 0; n < nSlices-blockSize; n++) zPosArray[n] = n*zStep;
        this.zStep = zStep;
    }

    public float computeZcorrelation(ImageStack ims){
        float zPos = 0;
        double[] zCorrArray = new double[ims.getSize()-blockSize];
        float[] zCorrArrayFloat = new float[ims.getSize()-blockSize];
//        IJ.log("Size of target block: "+ims.getSize());
//        IJ.log("Array size: "+zCorrArray.length);


        for (int z0=0; z0<(ims.getSize()-blockSize); z0++){
            zCorrArray[z0] = getZcorrelation(imsRef, ims.crop(0,0, z0, width, height, blockSize));
            zCorrArrayFloat[z0] = (float) zCorrArray[z0];
//            IJ.log("Value: "+zCorrArray[z0]);
        }

        float cropLevel = 0.7f;
        fitting = new Fit1DGaussian(normalizeArray(zCorrArrayFloat));
        fitting.cropDataArray(cropLevel);
        float[] fitResults = fitting.calculate();
        double[][] modelArray = fitting.fittedCurve();
        defocusPlot.setColor(Color.black);
        defocusPlot.add("line", zPosArray, modelArray[0]);
        defocusPlot.setColor(Color.red);
        defocusPlot.add("line", zPosArray, modelArray[1]);
//        defocusPlot.add("line", zPosArray, normalizeArray(zCorrArray));
        defocusPlot.show();

//        IJ.log("--- Initial guesses results ---");
//        IJ.log("Initial x0: "+(fitting.initX0)*zStep+" nm");
//        IJ.log("Initial Sigma: "+(fitting.initSigma)*zStep+" nm");

        IJ.log("--- Fit results ---");
        IJ.log("Amp: "+fitResults[0]);
        IJ.log("x0: "+(fitResults[1])*zStep+" nm");
        IJ.log("Sigma: "+(fitResults[2])*zStep+" nm");
        IJ.log("BG: "+fitResults[3]);

        return (fitResults[1])*zStep;
    }


    public double getZcorrelation(ImageStack imsRef, ImageStack ims){

        int nPixels = imsRef.getWidth()*imsRef.getHeight()*imsRef.getSize();
        float[] arrayRef = new float[nPixels];
        float[] array = new float[nPixels];
        arrayRef = imsRef.getVoxels(0, 0, 0, imsRef.getWidth(), imsRef.getHeight(), imsRef.getSize(), arrayRef);
        array = ims.getVoxels(0, 0, 0, imsRef.getWidth(), imsRef.getHeight(), imsRef.getSize(), array);

        float[] arrayProduct = multiply(arrayRef, array);
        return getSumValue(arrayProduct)/nPixels;

    }

    public double[] normalizeArray(double[] array){
        double max = getMaxValue(array)[1];
        double min = getMinValue(array)[1];

        double[] normArray = new double[array.length];
        for (int i=0; i<array.length; i++){
            normArray[i] = (array[i]-min)/(max-min);
        }

        return normArray;
    }

    public float[] normalizeArray(float[] array){
        float max = getMaxValue(array)[1];
        float min = getMinValue(array)[1];

        float[] normArray = new float[array.length];
        for (int i=0; i<array.length; i++){
            normArray[i] = (array[i]-min)/(max-min);
        }

        return normArray;
    }


}
