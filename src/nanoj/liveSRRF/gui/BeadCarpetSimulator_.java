package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import nanoj.core.java.io.LoadNanoJTable;
import nanoj.core2.NanoJPrefs;

import static nanoj.core.java.imagej.ResultsTableTools.dataMapToResultsTable;
import static nanoj.core2.NanoJFHT.*;


import java.io.IOException;
import java.util.Map;
import java.util.Random;

public class BeadCarpetSimulator_ implements PlugIn {

    public void run(String arg) {

//        String[] convMethods = new String[2];
//        int n = 0;
//        convMethods[n++] = "ImageJ Convolver class";
//        convMethods[n++] = "NanoJ FHT";

        NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
        String[] simulationPattern = new String[]{"random pattern", "grid pattern"};

        GenericDialog gd = new GenericDialog("Bead carpet simulator");
        gd.addChoice("Simulated pattern: ", simulationPattern, prefs.get("simPattern", simulationPattern[0]));
        gd.addNumericField("Image size (in nm)", prefs.get("imageSize", 5000), 1);
        gd.addNumericField("Bead diameter (in nm)", prefs.get("beadDiameter", 40), 1);
        gd.addNumericField("Bead density (in #/um^2)", prefs.get("beadDensity", 1.2f), 2);
        gd.addNumericField("Pixel size (in nm)", prefs.get("simPixelSize", 4), 1);
        gd.addCheckbox("Simulate displacement", prefs.get("doDisplacement", false));
        gd.addNumericField("Pixel size in displacement map: ", prefs.get("pixelSizeDisplacement", 50),2);
//        gd.addChoice("Convolution method", convMethods, convMethods[0]); // TODO: implement FHT convolution
        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        String simPattern = gd.getNextChoice();
        float imageSize = (float) gd.getNextNumber(); // nm
        float beadDiameter = (float) gd.getNextNumber(); // nm
        float beadDensity = (float) gd.getNextNumber(); // beads / um^2
        float simPixelSize = (float) gd.getNextNumber(); // nm
        boolean doDisplacement = gd.getNextBoolean();
        float pixelSizeDisplacement = (float) gd.getNextNumber();
//        String convMethodChosen = gd.getNextChoice();

        prefs.set("simPattern", simPattern);
        prefs.set("imageSize", imageSize);
        prefs.set("beadDiameter", beadDiameter);
        prefs.set("beadDensity", beadDensity);
        prefs.set("simPixelSize", simPixelSize);
        prefs.set("doDisplacement", doDisplacement);
        prefs.set("pixelSizeDisplacement", pixelSizeDisplacement);

        int nBeads = (int) (imageSize/1000 * imageSize/1000 * beadDensity);
        if (simPattern.equals("grid pattern")){
            nBeads = (int) Math.pow(Math.ceil(Math.sqrt(nBeads)),2);
        }
        IJ.log("\\Clear");
        IJ.log("------------------");
        IJ.log("Number of beads: "+nBeads);
        IJ.log("Bead diameter (nm): "+beadDiameter);

        Random random = new Random();
        int w = (int) (imageSize/simPixelSize);
        ImageStack imsDispXY;
        FloatProcessor fpUX = null;
        FloatProcessor fpUY = null;

        if (doDisplacement){
            imsDispXY = getDisplacementMapFromNanoJtable();
            if (imsDispXY == null) return;

            fpUX = imsDispXY.getProcessor(1).convertToFloatProcessor();
            fpUY = imsDispXY.getProcessor(2).convertToFloatProcessor();
            fpUX.setInterpolationMethod(fpUX.BICUBIC);
            fpUY.setInterpolationMethod(fpUY.BICUBIC);

            ImagePlus impDispXY = new ImagePlus("Displacement map (in pixels)", imsDispXY);
            impDispXY.show();
        }

        float[][] xyCoord = new float[nBeads][2]; // TODO: export the true localizations?
        float[][] xyCoordDisplaced = new float[nBeads][2]; // TODO: export the true localizations?

        float[][] pixels = new float[w][w];
        float[][] pixelsDisplaced = new float[w][w];

        IJ.log("Generating coordinates...");
        int x,y;
        for (int i = 0; i < nBeads; i++) {

            if (simPattern.equals("random pattern")) {
                xyCoord[i][0] = imageSize * random.nextFloat(); // in nm
                xyCoord[i][1] = imageSize * random.nextFloat(); // TODO: check distances between beads?
            }
            else {
                float beadRow = (float) Math.sqrt(nBeads);
                float xi = (int) (i / beadRow);
                float yi = (i - xi*beadRow);
                xyCoord[i][0] = imageSize * (xi + 0.5f) / beadRow; // in nm
                xyCoord[i][1] = imageSize * (yi + 0.5f) / beadRow; // TODO: check distances between beads?
            }

            x = Math.round(Math.max(Math.min(w*xyCoord[i][0]/imageSize, w-1),0));
            y = Math.round(Math.max(Math.min(w*xyCoord[i][1]/imageSize, w-1),0));

            pixels[x][y] = 1;
            if (doDisplacement) {
                xyCoordDisplaced[i][0] = xyCoord[i][0] + pixelSizeDisplacement * (float) fpUX.getInterpolatedPixel(fpUX.getWidth()*x/w, fpUX.getHeight()*y/w);
                xyCoordDisplaced[i][1] = xyCoord[i][1] + pixelSizeDisplacement * (float) fpUY.getInterpolatedPixel(fpUX.getWidth()*x/w, fpUX.getHeight()*y/w); // TODO: check distances between beads?
                x = Math.round(Math.max(Math.min(w*xyCoordDisplaced[i][0]/imageSize, w-1),0));
                y = Math.round(Math.max(Math.min(w*xyCoordDisplaced[i][1]/imageSize, w-1),0));
                pixelsDisplaced[x][y] = 1;
            }
        }

        ImageStack ims = new ImageStack(w,w);
        FloatProcessor fp = new FloatProcessor(pixels);
        ims.addSlice(fp);
        ImagePlus imp = new ImagePlus("Ground truth bead carpet", ims);
        IJ.run(imp, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
        imp.show();
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");

        FloatProcessor fpDisplaced = null;
        if (doDisplacement) {
            ImageStack imsDisplaced = new ImageStack(w, w);
            fpDisplaced = new FloatProcessor(pixelsDisplaced);
            imsDisplaced.addSlice(fpDisplaced);
            ImagePlus impDisplaced = new ImagePlus("Ground truth bead carpet - displaced", imsDisplaced);
            IJ.run(impDisplaced, "Set Scale...", "distance=1 known=" + simPixelSize + " unit=nm");
            impDisplaced.show();
            IJ.run(impDisplaced, "Enhance Contrast", "saturated=0.35");
        }

        // Generating the bead kernel
        IJ.log("Generating bead kernel...");
        float[] beadKernel = getBeadKernel(beadDiameter, simPixelSize);
        int kernelSize = (int) Math.sqrt(beadKernel.length);
        IJ.log("Kernel size: "+kernelSize);
        ImageStack imsKernel = new ImageStack(kernelSize, kernelSize);
        imsKernel.addSlice(new FloatProcessor(kernelSize, kernelSize, beadKernel));
        ImagePlus impKernel = new ImagePlus("Bead kernel", imsKernel);
        IJ.run(impKernel, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
        impKernel.show();
        IJ.run(impKernel, "Enhance Contrast", "saturated=0.35");

        // Convolving
        FloatProcessor fpConv = fp.duplicate().convertToFloatProcessor();
        IJ.log("Convolving...");
        Convolver cv = new Convolver();
        cv.setNormalize(false);
        cv.convolve(fpConv, beadKernel, kernelSize, kernelSize);

        ImageStack imsConv = new ImageStack(w,w);
        imsConv.addSlice(fpConv);
        ImagePlus impConv = new ImagePlus("Simulated bead carpet", imsConv);
        IJ.run(impConv, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
        impConv.show();
        IJ.run(impConv, "Enhance Contrast", "saturated=0.35");

        if (doDisplacement){
            FloatProcessor fpConvDisplaced = fpDisplaced.duplicate().convertToFloatProcessor();
            cv.convolve(fpConvDisplaced, beadKernel, kernelSize, kernelSize);
            ImageStack imsConvDisplaced = new ImageStack(w,w);
            imsConvDisplaced.addSlice(fpConvDisplaced);
            ImagePlus impConvDisplaced = new ImagePlus("Simulated bead carpet - Displaced", imsConvDisplaced);
            IJ.run(impConvDisplaced, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
            impConvDisplaced.show();
            IJ.run(impConvDisplaced, "Enhance Contrast", "saturated=0.35");
        }


        IJ.log("------------------");
        IJ.log("All done.");


    }

    final static float[] getBeadKernel(float beadDiameter, float pixelSize){

        int kernelSize = (int) (beadDiameter/pixelSize) + 2;
//        IJ.log("Kernel size (in function): "+kernelSize);
//        IJ.log("Kernel size modulus: "+(kernelSize % 2));

        if ((kernelSize % 2) == 0) kernelSize++; // make it odd

//        IJ.log("Kernel size (in function): "+kernelSize);

        float[] kernel = new float[kernelSize*kernelSize];
        int xy0 = (kernelSize-1)/2;

        for (int i = 0; i < kernelSize; i++) {
            float xi = pixelSize*(i - xy0);
            for (int j = 0; j < kernelSize; j++) {
                float yi =  pixelSize*(j - xy0);
                float r = (float) Math.sqrt(xi*xi + yi*yi);
                if (r < 0.5*beadDiameter){
                    kernel[i*kernelSize + j] = (float) Math.sqrt(0.25*beadDiameter*beadDiameter - r*r);}
                else
                    kernel[i*kernelSize + j] = 0;
            }
        }

        return kernel;
    }



    final static ImageStack getDisplacementMapFromNanoJtable(){

        // ---- Getting calibration data from the NanoJ table ----
        IJ.log("Getting displacement map...");
        String displacementTablePath = IJ.getFilePath("Choose displacement table to load...");
        if (displacementTablePath == null) return null;

        double[] uxDis, uyDis;
        uxDis = null;
        uyDis = null;

        Map<String, double[]> displacementTable;
        try {
            displacementTable = new LoadNanoJTable(displacementTablePath).getData(); // TODO: accept other data format??
            uxDis = displacementTable.get("ux1");
            uyDis = displacementTable.get("uy1");
            ResultsTable rt = dataMapToResultsTable(displacementTable);
            rt.show("Displacement-Table");
        } catch (IOException e) {
            IJ.log("Catching exception...");
            e.printStackTrace();
        }

        int dispMapWidth = (int) Math.sqrt(uxDis.length);

        double[] dispMag = new double[uxDis.length];
        double[] dispAngle = new double[uxDis.length];
        for (int i = 0; i < uxDis.length; i++) {
            dispMag[i] = Math.sqrt(uxDis[i]*uxDis[i] + uyDis[i]*uyDis[i]);
            dispAngle[i] = Math.atan(uyDis[i]/uxDis[i]);
        }

        FloatProcessor fpUX = new FloatProcessor(dispMapWidth, dispMapWidth, uxDis);
        FloatProcessor fpUY = new FloatProcessor(dispMapWidth, dispMapWidth, uyDis);
        FloatProcessor fpMagDisp = new FloatProcessor(dispMapWidth, dispMapWidth, dispMag);
        FloatProcessor fpAngleDisp = new FloatProcessor(dispMapWidth, dispMapWidth, dispAngle);

        ImageStack imsDispXY = new ImageStack(dispMapWidth, dispMapWidth);
        imsDispXY.addSlice(fpUX);
        imsDispXY.addSlice(fpUY);
        imsDispXY.addSlice(fpMagDisp);
        imsDispXY.addSlice(fpAngleDisp);

        return imsDispXY;
    }


}
