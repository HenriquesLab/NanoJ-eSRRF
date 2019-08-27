package nanoj.liveSRRF.gui;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;

public class JavaTest_ implements PlugIn {
    public void run(String arg) {

        // Get raw data
        ImagePlus imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        ImageStack ims = imp.getImageStack();
        ImageStack imsReshaped = reshapeImageStack3Dv3(ims);
        ImagePlus impReshaped = new ImagePlus("Reshaped", imsReshaped);
        impReshaped.show();

        FloatProcessor fp = imsReshaped.getProcessor(1).convertToFloatProcessor();
        fp.add(6);



        IJ.log("All done.");

    }

    // --- Reshape the data for 3D-SRRF ---
    public ImageStack reshapeImageStack3D(ImageStack imsIn) {


        ImageStack ims = imsIn.duplicate();
        int widthS = ims.getWidth() / 3;
        int heightS = ims.getHeight() / 3;
        int nPlanesThreeD = 9;
        ImageStack reshapedIms = new ImageStack(widthS * heightS, nPlanesThreeD, ims.getSize());

        float[] pixels = new float[widthS * heightS];
        int x, y;
        for (int f = 0; f < ims.getSize(); f++) {
            for (int z = 0; z < nPlanesThreeD; z++) {
                y = (z / 3) * heightS;
                x = (z - 3 * (z / 3)) * heightS;
                IJ.log("z " + z);
                IJ.log("x/y: " + x + "/" + y);
                pixels = ims.getVoxels(x, y, f, widthS, heightS, 1, pixels);
                reshapedIms.setVoxels(0, z, f, widthS * heightS, 1, 1, pixels);
            }

        }
        return reshapedIms;
    }


    // --- Reshape the data for 3D-SRRF ---
    public ImageStack reshapeImageStack3Dv2(ImageStack imsTemp) {

//        ImageStack imsTemp = ims.duplicate();
        int width = imsTemp.getWidth();
        int height = imsTemp.getHeight();
        int widthS = width/3;
        int heightS = height/3;
        int nPlanesThreeD = 9;
        ImageStack reshapedIms = new ImageStack(widthS*heightS, nPlanesThreeD, imsTemp.getSize());

        float[] pixels;
        float[] pixelsReshaped = new float[widthS*heightS*nPlanesThreeD];

        int z, xL, yL;
        for (int f = 0; f < imsTemp.getSize(); f++) {
            IJ.log("Frame "+f);
            pixels = (float[]) imsTemp.getProcessor(f+1).getPixels();
            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    z = (y/heightS)*3 + (x/widthS);
                    yL = y - (z/3)*heightS;
                    xL = x - (z - 3*(z/3))*widthS;
//                    IJ.log("z " + z);
//                    IJ.log("x/y: " + xL + "/" + yL);
                    pixelsReshaped[xL + yL*widthS + z*widthS*heightS] = pixels[x + y*width];
//                    pixelsReshaped[xL + yL*widthS + z*widthS*heightS] = (float) ims.duplicate().getVoxel(x, y, f);

                }
            }
            FloatProcessor fp = new FloatProcessor(widthS*heightS, nPlanesThreeD, pixelsReshaped);
            reshapedIms.setProcessor(fp, f+1);
        }

        return reshapedIms;

    }

    // --- Reshape the data for 3D-SRRF ---
    public ImageStack reshapeImageStack3Dv3(ImageStack imsTemp) {

//        ImageStack imsTemp = ims.duplicate();
        int width = imsTemp.getWidth();
        int height = imsTemp.getHeight();
        int widthS = width/3;
        int heightS = height/3;
        int nPlanesThreeD = 9;
        ImageStack reshapedIms = new ImageStack(widthS*heightS, nPlanesThreeD, imsTemp.getSize());


        int z, xL, yL, xG, yG;
        for (int f = 0; f < imsTemp.getSize(); f++) {
//            IJ.log("Frame "+f);
            FloatProcessor fp = imsTemp.getProcessor(f+1).convertToFloatProcessor();
            float[] pixelsReshaped = new float[widthS*heightS*nPlanesThreeD];

            for (int i = 0; i < widthS*heightS*nPlanesThreeD; i++) {
                z = i/(widthS*heightS);
                yL = (i - z*widthS*heightS)/widthS; // local coordinates
                xL = i - z*widthS*heightS - yL*widthS;
                xG = xL + (z - 3*(z/3))*widthS;
                yG = yL + (z/3)*heightS;
//                pixelsReshaped[i] = (float) imsTemp.getVoxel(xG, yG, f);
//                pixelsReshaped[i] = pixels[xG + yG*width];
                pixelsReshaped[i] = fp.getf(xG, yG);
            }
            FloatProcessor fpOut = new FloatProcessor(widthS*heightS, nPlanesThreeD, pixelsReshaped);
            reshapedIms.setProcessor(fpOut, f+1);
        }

        return reshapedIms;

    }
}
