package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

import static java.lang.Math.pow;

public class CalculateGainAndOffset_  implements PlugIn {

    @Override
    public void run(String s) {

        String dirPath = IJ.getDirectory("Choose directory to open...");
        if (dirPath == null) return;

        // need to add extension question here

        ArrayList<String> filesToProcess = this.getFilesToProcess(dirPath, ".tif");

        ImageStack imsMean = null;
        ImageStack imsVar = null;

        for (String fileName: filesToProcess) {

            IJ.log(fileName);
            IJ.run("Bio-Formats Windowless Importer", "open=["+fileName+"]");
            ImagePlus imp = IJ.getImage();
//            ImagePlus imp = IJ.openImage(fileName);

            ImageStack ims = imp.getImageStack();

            int w = ims.getWidth();
            int h = ims.getHeight();

            if (imsMean == null) {
                imsMean = new ImageStack(w, h);
                imsVar = new ImageStack(w, h);
            }
            else {
                assert (imsMean.getWidth() == w);
                assert (imsMean.getHeight() == h);
            }

            FloatProcessor[] fps = this.calculateMeanAndVarianceProjection(ims);
            imsMean.addSlice(fps[0]);
            imsVar.addSlice(fps[1]);
        }

        new ImagePlus("Mean", imsMean).show();
        new ImagePlus("Var", imsVar).show();
    }

    private ArrayList<String> getFilesToProcess(String dirPath, String extension) {

        File directory = new File(dirPath);
        if(!directory.exists()){
            return null;
        }

        ArrayList<String> filesToProcess = new ArrayList<String>();

        File[] listOfFilesInFolder = directory.listFiles();
        Arrays.sort(listOfFilesInFolder);

        for (int i=0;i<listOfFilesInFolder.length;i++){
            IJ.showProgress(i+1, listOfFilesInFolder.length);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return null;
            }

            if(listOfFilesInFolder[i].isFile() && listOfFilesInFolder[i].toString().endsWith(extension)) {
                String fPath = listOfFilesInFolder[i].getPath();

                filesToProcess.add(fPath);
            }
        }

        return filesToProcess;
    }

    private FloatProcessor[] calculateMeanAndVarianceProjection(ImageStack ims) {
        int w = ims.getWidth();
        int h = ims.getHeight();

        float[] mean = new float[w*h];
        float[] var  = new float[w*h];

        for (int s=1; s<=ims.getSize(); s++) {
            FloatProcessor frame = ims.getProcessor(s).convertToFloatProcessor();
            float[] pixels = (float[]) frame.getPixels();

            for (int n=0; n<pixels.length; n++) {
                mean[n] += (pixels[n]-mean[n]) / s;
                var[n] += pow(pixels[n]-mean[n], 2) / s; // check on this
            }
        }

        FloatProcessor fpMean = new FloatProcessor(w, h, mean);
        FloatProcessor fpVar = new FloatProcessor(w, h, var);

        return new FloatProcessor[] {fpMean, fpVar};
    }
}
