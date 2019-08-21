package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

public class LinearRescaling_ implements PlugIn {

    private int width;
    private int height;

    public void run(String arg) {
        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        IJ.run(imp,"32-bit", "");
        imp.show();


        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("y = a + b*x");

        int nSlices = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();
        String imageTitle = imp.getTitle();

        String label;

        ImageStack imsInput = imp.getImageStack();
        ImageStack imsOutput = new ImageStack(width, height, nSlices);

        ImageProcessor ipRef = imsInput.getProcessor(1);
        //FloatProcessor fpRef = ipRef.convertToFloatProcessor();
        float[] pixelsRef = (float[]) ipRef.getPixels();

        imsOutput.setProcessor(ipRef, 1);
        imsOutput.setSliceLabel(imsInput.getSliceLabel(1), 1);

        for (int s = 2; s <= nSlices; s++) {
            ImageProcessor ip = imsInput.getProcessor(s);

            float[] pixels = (float[]) ip.getPixels();
            float[] ab = getLinearRegressionParameters(pixelsRef, pixels);
            IJ.log("a="+ab[0]+" / b="+ab[1]);

            imsOutput.setProcessor(getRescaledFP(ab, pixels), s);
            imsOutput.setSliceLabel(imsInput.getSliceLabel(s), s);

        }

        ImagePlus impOutput = new ImagePlus(imageTitle + " - Linearly rescaled", imsOutput);
        impOutput.show();
    }


    public static float[] getLinearRegressionParameters(float[] floatArrayRef, float[] floatArray){

        int n = floatArrayRef.length;
        float x_bar = 0;
        float y_bar = 0;
        float x2 = 0;
        float xy = 0;

        for (int i=0; i<n; i++){
            x_bar += floatArray[i]/(float) n;
            y_bar += floatArrayRef[i]/(float) n;
            x2 += floatArray[i]*floatArray[i];
            xy += floatArray[i]*floatArrayRef[i];
        }

        float[] ab = new float[2];
        ab[0] = (y_bar*x2 - x_bar*xy)/(x2 - n*x_bar*x_bar);
        ab[1] = (xy - n*x_bar*y_bar)/(x2-n*x_bar*x_bar);
        return ab;
    }


    private FloatProcessor getRescaledFP(float[] ab, float[] floatArray){

        int n = floatArray.length;
        float[] floatArrayOut = new float[n];
        for (int i=0; i<n; i++){
            floatArrayOut[i] = ab[1]*floatArray[i] + ab[0];
        }
        return new FloatProcessor(width, height, floatArrayOut);
    }


}
