package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.plugin.PlugIn;
import ij.plugin.filter.Convolver;
import ij.process.FloatProcessor;
import static nanoj.core2.NanoJFHT.*;


import java.util.Random;

public class BeadCarpetSimulator_ implements PlugIn {

    public void run(String arg) {

//        String[] convMethods = new String[2];
//        int n = 0;
//        convMethods[n++] = "ImageJ Convolver class";
//        convMethods[n++] = "NanoJ FHT";

        GenericDialog gd = new GenericDialog("Bead carpet simulator");
        gd.addNumericField("Image size (in nm)", 5000, 1);
        gd.addNumericField("Bead diameter (in nm)", 40, 1);
        gd.addNumericField("Bead density (in #/um^2)", 1.2, 2);
        gd.addNumericField("Pixel size (in nm)", 4, 1);
//        gd.addChoice("Convolution method", convMethods, convMethods[0]);

        gd.showDialog();

        if (gd.wasCanceled()) {
            return;
        }

        float imageSize = (float) gd.getNextNumber(); // nm
        float beadDiameter = (float) gd.getNextNumber(); // nm
        float beadDensity = (float) gd.getNextNumber(); // beads / um^2
        float simPixelSize = (float) gd.getNextNumber(); // nm
//        String convMethodChosen = gd.getNextChoice();

        int nBeads = (int) (imageSize/1000 * imageSize/1000 * beadDensity);
        IJ.log("\\Clear");
        IJ.log("------------------");
        IJ.log("Number of beads: "+nBeads);
        IJ.log("Bead diameter (nm): "+beadDiameter);

        Random random = new Random();
        int w = (int) (imageSize/simPixelSize);

        float[][] xyCoord = new float[nBeads][2];
        float[][] pixels = new float[w][w];

        IJ.log("Generating coordinates...");
        for (int i = 0; i < nBeads; i++) {
            xyCoord[i][0] = imageSize*random.nextFloat(); // in nm
            xyCoord[i][1] = imageSize*random.nextFloat(); // TODO: check distances between beads?
            pixels[Math.round(w*xyCoord[i][0]/imageSize)][Math.round(w*xyCoord[i][1]/imageSize)] = 1;
        }

        ImageStack ims = new ImageStack(w,w);
        ims.addSlice(new FloatProcessor(pixels));
        ImagePlus imp = new ImagePlus("Ground truth bead carpet", ims);
        IJ.run(imp, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
        imp.show();
        IJ.run(imp, "Enhance Contrast", "saturated=0.35");

        // Generating the bead kernel
        IJ.log("Generating bead kernel...");
        FloatProcessor fp = new FloatProcessor(pixels);
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
        FloatProcessor fpConv;
//        if (convMethodChosen.equals("ImageJ Convolver class")){
            fpConv = fp.duplicate().convertToFloatProcessor();
        IJ.log("Convolving...");
        Convolver cv = new Convolver();
        cv.setNormalize(false);
        cv.convolve(fpConv, beadKernel, kernelSize, kernelSize);
//    }
//        else {
//            fpConv = fhtConvolution(fp, new FloatProcessor(kernelSize, kernelSize, beadKernel));
//        }


        ImageStack imsConv = new ImageStack(w,w);
        imsConv.addSlice(fpConv);
        ImagePlus impConv = new ImagePlus("Simulated bead carpet", imsConv);
        IJ.run(impConv, "Set Scale...", "distance=1 known="+simPixelSize+" unit=nm");
        impConv.show();
        IJ.run(impConv, "Enhance Contrast", "saturated=0.35");
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

//    final static FloatProcessor fhtConvolution(FloatProcessor fp, FloatProcessor fpKernel){ // TODO: not done yet
//        FloatProcessor fpFHT = forwardFHT(fp);
//        FloatProcessor fpKernelFHT = forwardFHT(fpKernel);
//
//        return fpFHT;
//    }
//
//    final static FloatProcessor resizeFloatProcessorToNextPow2(FloatProcessor fp){
//
//        int w = fp.getWidth();
//        int h = fp.getHeight();
//
//        int bitW = (int) Math.ceil(Math.log(w)/Math.log(2));
//        if (Math.pow(2,bitW) == w) bitW++;
//
//        int bitH = (int) Math.ceil(Math.log(h)/Math.log(2));
//        if (Math.pow(2,bitH) == h) bitH++;
//
//
//
//
//        return fpFHT;
//    }

}
