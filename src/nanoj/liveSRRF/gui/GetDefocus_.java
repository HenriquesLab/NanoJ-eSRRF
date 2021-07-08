package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Plot;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.EntropyMeasure;



public class GetDefocus_ implements PlugIn {
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());


    public void run(String s) {

        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        imp.show();

        FileInfo imageInfo;
        imageInfo = imp.getFileInfo();
        double zStep = 1000*imageInfo.pixelDepth; //nm
        zStep = (double) Math.round(zStep*100d)/100d;

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Get defocus from Image splitter:");
        gd.addNumericField("Number of split/axis:", prefs.get("nSplits", 3), 0);
        gd.addNumericField("Crop level:", prefs.get("cropLevel", 0.5f), 2);

        gd.showDialog();
        if (gd.wasCanceled()) return;

        int nImageSplits = (int) gd.getNextNumber();
        float cropLevel = (float) gd.getNextNumber();

        prefs.set("nSplits", (float) nImageSplits);
        prefs.set("cropLevel", (float) cropLevel);

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------------------------------------------------------");
        IJ.log("z-step: "+zStep+" nm");
        IJ.log("Level for cropping: "+cropLevel);

        // ---- Getting input data ----
        ImageStack ims = imp.duplicate().getImageStack().convertToFloat();

        int width = ims.getWidth();
        int height = ims.getHeight();
        int nSlices = ims.getSize();

        int borderCrop = 15;

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