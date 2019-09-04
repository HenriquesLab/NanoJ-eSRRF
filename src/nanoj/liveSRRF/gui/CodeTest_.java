package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core.java.array.ArrayMath;
import nanoj.liveSRRF.EntropyMeasure;
import nanoj.liveSRRF.ShannonEntropyMeasure;

import static nanoj.liveSRRF.StackProjections.calculateAverage;
import static nanoj.liveSRRF.StackProjections.calculateMIP;
import static nanoj.liveSRRF.gui.GetSpatialCalibrationMFMdata_.getMasksFromIntensityPercentile;

public class CodeTest_ implements PlugIn {
    public void run(String arg) {

        // ---- Getting input data ----
        ImagePlus imp = WindowManager.getCurrentImage();
        ImageStack ims = imp.duplicate().getImageStack().convertToFloat();

        int width = ims.getWidth();
        int height = ims.getHeight();
        int nSlices = ims.getSize();
        int nImageSplits = 3;
        int borderCrop = 15;

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000*imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep*100d)/100d;

        double[] zPos = new double[nSlices];
        for (int n = 0; n < nSlices; n++) zPos[n] = n*zStep;

        // ---- Initialize the variables and the Log ----
        int cropSizeX = width/nImageSplits - 2*borderCrop; // int division rounds down automatically, ignoring the remainder
        int cropSizeY = height/nImageSplits - 2*borderCrop;

        ImageStack[] imsArray = new ImageStack[nImageSplits*nImageSplits];

        // ---- Cropping loop based on nominal xyz positions ----
        IJ.log("------------");
        IJ.log("Cropping FOV...");
        int x,y;
        for (int j = 0; j < nImageSplits; j++){
            for (int i = 0; i < nImageSplits; i++){
                x = Math.round((float) width / nImageSplits * i) + borderCrop;
                y = Math.round((float) height / nImageSplits * j) + borderCrop;
//                IJ.log("x/y/z: "+x+"/"+y+"/"+z);
                imsArray[i + j*nImageSplits] = ims.duplicate().crop(x, y, 0, cropSizeX, cropSizeY, nSlices);
            }
        }
        Plot defocusPlot = new Plot("Defocus plot", "Z (nm)", "Entropy (AU)");

        for (int n = 0; n < nImageSplits*nImageSplits; n++) {
            IJ.log("Slice #"+n);
//            ShannonEntropyMeasure entropyMeasure = new ShannonEntropyMeasure(imsArray[n]);
//            double[] entropyArray = entropyMeasure.EntropyArray;
            EntropyMeasure entropyMeasure = new EntropyMeasure(imsArray[n]);
            double[] thisEntropyArray = entropyMeasure.entropyArray;
            defocusPlot.add("line", zPos, thisEntropyArray);
        }

        IJ.log("THis try "+Math.log(0));
        defocusPlot.show();
        IJ.log("All done");

    }


}
