package nanoj.liveSRRF;

import com.sun.jna.platform.win32.WinBase;
import ij.ImageStack;
import ij.gui.Plot;
import ij.process.FloatProcessor;
import nanoj.core.java.image.analysis.FRC;

public class FHT3D_Utilities {

    public static ImageStack reshape3DFHTforInterpolation(ImageStack ims, int magnification){

        int nSlices = ims.getSize();
        ImageStack imsOut = new ImageStack(ims.getWidth()*magnification, ims.getHeight()*magnification);
        FloatProcessor fpInt;
        for (int z = 0; z < nSlices*magnification; z++) {
            if (z < nSlices/2){
//                IJ.log("z position: "+z+" ---> mode 1");
                FloatProcessor fp = ims.getProcessor(z+1).convertToFloatProcessor();
                fpInt = fhtSpaceInterpolation(fp, magnification);
            }
            else if (z < (nSlices*magnification - nSlices/2)){
//                IJ.log("z position: "+z+" ---> mode 2");
                float[] pixels = new float[ims.getWidth()*magnification * ims.getHeight()*magnification];
                fpInt = new FloatProcessor(ims.getWidth()*magnification, ims.getHeight()*magnification, pixels);
            }
//            else if (z >= (nSlices*magnification - nSlices/2)){
            else {
//                IJ.log("z position: "+z+" ---> mode 3");
                int zOffset = z - nSlices*(magnification - 1);
//                IJ.log("z offset: "+zOffset);
                FloatProcessor fp = ims.getProcessor(zOffset+1).convertToFloatProcessor();
                fpInt = fhtSpaceInterpolation(fp, magnification);
            }

            imsOut.addSlice(fpInt);
        }

        return imsOut;

    }


    public static FloatProcessor fhtSpaceInterpolation(FloatProcessor fpFHT, int intFactor){

        int w = fpFHT.getWidth();
        int h = fpFHT.getHeight();
        int wInt = w*intFactor;
        int hInt = h*intFactor;
        float[] pixelsFHT = (float[]) fpFHT.duplicate().getPixels();
        float[] pixelsFHTint = new float[wInt * hInt];
        float intF2 = (float) intFactor*intFactor;

        for (int p = 0; p < wInt; p++) {
            for (int q = 0; q < hInt; q++) {
                if (p <= (w / 2)) {
                    if (q <= (h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[p * w + q];
                    if ((q > (h / 2)) && (q <= (hInt - h / 2))) pixelsFHTint[p * wInt + q] = 0;
                    if (q > (hInt - h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[p*w + q + h - hInt];
                }

                if ((p > (w / 2)) && (p <= wInt - w / 2) && (q <= (h / 2))) pixelsFHTint[p * wInt + q] = 0;

                if (p > (wInt - w / 2)) {
                    if (q <= (h / 2)) pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[(p + w - wInt) * w + q];
                    if ((q > (h / 2)) && (q <= (hInt - h / 2))) pixelsFHTint[p * wInt + q] = 0;
                    if (q > (hInt - h / 2))
                        pixelsFHTint[p * wInt + q] = intF2 * pixelsFHT[(p + w - wInt) * w + q + h - hInt];
                }
            }
        }

        return new FloatProcessor(wInt, hInt, pixelsFHTint);
    }

    public static FloatProcessor mirrorPadding(FloatProcessor fp){

        int w = fp.getWidth();
        int h = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();
        float[] pixels3x3 = new float[3*w*3*h];

        int bitW = (int) Math.ceil(Math.log(w)/Math.log(2));
        if (Math.pow(2,bitW) == w) bitW++;
        bitW = (int) Math.pow(2, bitW);

        int bitH = (int) Math.ceil(Math.log(h)/Math.log(2));
        if (Math.pow(2,bitH) == h) bitH++;
        bitH = (int) Math.pow(2, bitH);

        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                pixels3x3[j*3*w + i] = pixels[(w-1-j)*w + (h-1-i)];
                pixels3x3[(j+h)*3*w + i] = pixels[j*w + (h-1-i)];
                pixels3x3[(j+2*h)*3*w + i] = pixels[(w-1-j)*w + (h-1-i)];

                pixels3x3[j*3*w + i+h] = pixels[(w-1-j)*w + i];
                pixels3x3[(j+h)*3*w + i+h] = pixels[j*w + i];
                pixels3x3[(j+2*h)*3*w + i+h] = pixels[(w-1-j)*w + i];

                pixels3x3[j*3*w + i+2*h] = pixels[(w-1-j)*w + (h-1-i)];
                pixels3x3[(j+h)*3*w + i+2*h] = pixels[j*w + (h-1-i)];
                pixels3x3[(j+2*h)*3*w + i+2*h] = pixels[(w-1-j)*w + (h-1-i)];
            }
        }

        FloatProcessor fpPadded = new FloatProcessor(3*w, 3*h, pixels3x3);

        int xROI = (3*w-bitW)/2; // TODO: currently works only if original image is even sized, check
        int yROI = (3*h-bitH)/2;

        fpPadded.setRoi(xROI, yROI, bitW, bitH);
        FloatProcessor fpCropped = fpPadded.crop().convertToFloatProcessor();

        return fpCropped;

    }

    public static ImageStack swapQuadrants3D(ImageStack ims){

        ImageStack imsSwapped = new ImageStack(ims.getWidth(), ims.getHeight());

        // Initialise the ImageStack
        for (int z = 0; z < ims.getSize(); z++) {
            float[] pixels = new float[ims.getWidth() * ims.getHeight()];
            FloatProcessor fp = new FloatProcessor(ims.getWidth(), ims.getHeight(), pixels);
            imsSwapped.addSlice(fp);
        }

        int blockSize = ims.getWidth()/2; // assumes that the ims is square so w=h=d
        int arrayLength = ims.getWidth()/2*ims.getHeight()/2*ims.getSize()/2;

        float[] pixelBlock = new float[arrayLength];
        // Swapping bottom block ------------
        pixelBlock = ims.getVoxels(0,0,0, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(blockSize, blockSize, blockSize, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(blockSize, 0,0, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(0, blockSize, blockSize, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(0, blockSize,0, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(blockSize, 0, blockSize, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(blockSize, blockSize,0, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(0, 0, blockSize, blockSize, blockSize, blockSize, pixelBlock);

        // Swapping top block ------------
//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(0,0,blockSize, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(blockSize, blockSize, 0, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(blockSize, 0,blockSize, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(0, blockSize, 0, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(0, blockSize,blockSize, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(blockSize, 0, 0, blockSize, blockSize, blockSize, pixelBlock);

//        pixelBlock = new float[arrayLength];
        pixelBlock = ims.getVoxels(blockSize, blockSize,blockSize, blockSize, blockSize, blockSize, pixelBlock);
        imsSwapped.setVoxels(0, 0, 0, blockSize, blockSize, blockSize, pixelBlock);

        return imsSwapped;
    }

    public static ImageStack[] getComplexFFT(ImageStack imsFHT){

        ImageStack[] imsFFTarray = new ImageStack[2];
        int imageSize = imsFHT.getSize();
        float[] pixelsReal;
        float[] pixelsImaginary;

        for (int z = 0; z < imageSize; z++) {
            pixelsReal = new float[imageSize * imageSize];
            pixelsImaginary = new float[imageSize * imageSize];
            for (int y = 0; y < imageSize; y++) {
                for (int x = 0; x < imageSize; x++) {
                    pixelsReal[x + imageSize*y] = 0.5f * (float) (imsFHT.getVoxel(x,y,z) + imsFHT.getVoxel(imageSize - x-1, imageSize - y-1, imageSize - z-1));
                    pixelsImaginary[x + imageSize*y] = 0.5f * (float) (imsFHT.getVoxel(x,y,z) - imsFHT.getVoxel(imageSize - x-1, imageSize - y-1, imageSize - z-1));
                }
            }
            imsFFTarray[0].addSlice(new FloatProcessor(imageSize, imageSize, pixelsReal));
            imsFFTarray[1].addSlice(new FloatProcessor(imageSize, imageSize, pixelsImaginary));
        }

        return imsFFTarray;
    }

    public static ImageStack getFFTmagnitude(ImageStack imsFHT){

        ImageStack imsFHTswapped = swapQuadrants3D(imsFHT);
        ImageStack imsFFTmagnitude = new ImageStack(imsFHT.getWidth(), imsFHT.getHeight());
        int imageSize = imsFHT.getSize();
        float[] pixels;

        for (int z = 0; z < imageSize; z++) {
            pixels = new float[imageSize * imageSize];
            for (int y = 0; y < imageSize; y++) {
                for (int x = 0; x < imageSize; x++) {
                    double h = imsFHTswapped.getVoxel(x,y,z);
                    double h_ = imsFHTswapped.getVoxel(imageSize - x-1, imageSize - y-1, imageSize - z-1);
                    pixels[x + imageSize*y] = (float) Math.sqrt(0.5*(h*h + h_*h_));
                }
            }
            imsFFTmagnitude.addSlice(new FloatProcessor(imageSize, imageSize, pixels));
        }

        return imsFFTmagnitude;
    }

    public static double getFSCresolution(ImageStack imsFHT1, ImageStack imsFHT2, float pixelSize){

        ImageStack imsFHT1swap = swapQuadrants3D(imsFHT1);
        ImageStack imsFHT2swap = swapQuadrants3D(imsFHT2);
        int imageSize = imsFHT1.getSize();

//        FRC.FRCCurveResult[] frcResults = new FRC.FRCCurveResult[imageSize/2+1]; // TODO: use the FRC Curve results

        double[] frcCurve = new double[imageSize/2];
        double[] frequencyArray = new double[imageSize/2];
        double[] threshold = new double[imageSize/2];
        double[] radiusArray = new double[imageSize/2];
        double[] mag1 = new double[imageSize/2];
        double[] mag2 = new double[imageSize/2];
        double[] numSumCount = new double[imageSize/2];

        int nAngles;
        double angleStep; // half circle
        double angle0;

        float R = 0.5f;
        int i = 0;
        while (R < imageSize/2.0f){
            float z = -R;
            while (z <= R){

                double r = Math.sqrt(R*R - z*z);
                if (r == 0) {
                    nAngles = 1;
                    angleStep = 0;
                    angle0 = 0;
                }
                else {
                    nAngles = (int) (Math.PI*r);
                    angleStep = Math.PI/nAngles; // half circle
                    angle0 = 0.5d/((int) r + 0.5d); // smallest angle aligned with a pixel, therefore this minimizes the amount of interpolation
//                    angle0 = 0;
                }

                FloatProcessor fp1 = imsFHT1swap.getProcessor((imageSize/2) + (int) (z+0.5f)).convertToFloatProcessor();
                FloatProcessor fp2 = imsFHT2swap.getProcessor((imageSize/2) + (int) (z+0.5f)).convertToFloatProcessor();
                FloatProcessor fp1_ = imsFHT1swap.getProcessor((imageSize/2) + (int) (-z+0.5f)).convertToFloatProcessor();
                FloatProcessor fp2_ = imsFHT2swap.getProcessor((imageSize/2) + (int) (-z+0.5f)).convertToFloatProcessor();

                fp1.setInterpolationMethod(fp1.BICUBIC);
                fp2.setInterpolationMethod(fp2.BICUBIC);
                fp1_.setInterpolationMethod(fp1_.BICUBIC);
                fp2_.setInterpolationMethod(fp2_.BICUBIC);

                for (int na = 0; na < nAngles; na++) {
                    double x = r*Math.cos(na*angleStep + angle0);
                    double y = r*Math.sin(na*angleStep + angle0);
                    double h1 = fp1.getInterpolatedPixel(x + imageSize/2.0d - 0.5, y + imageSize/2.0d - 0.5);
                    double h2 = fp2.getInterpolatedPixel(x + imageSize/2.0d - 0.5, y + imageSize/2.0d - 0.5);

                    x = r*Math.cos(Math.PI + na*angleStep + angle0);
                    y = r*Math.sin(Math.PI + na*angleStep + angle0);
                    double h1_ = fp1_.getInterpolatedPixel(x + imageSize/2.0d - 0.5, y + imageSize/2.0d - 0.5);
                    double h2_ = fp2_.getInterpolatedPixel(x + imageSize/2.0d - 0.5, y + imageSize/2.0d - 0.5);
                    frcCurve[i] += h1*h2 + h1_*h2_;
                    mag1[i] += h1*h1 + h1_*h1_;
                    mag2[i] += h2*h2 + h2_*h2_;
                    numSumCount[i] ++;
                }
                z += 1.0f;
            }

            radiusArray[i] = R;
//            frcResults[i] = new FRC.FRCCurveResult((int) (R+0.5), (int) numSumCount[i], frcCurve[i], mag1[i], mag2[i]);
            frcCurve[i] = frcCurve[i]/(Math.sqrt(mag1[i]*mag2[i]));
            threshold[i] = 1.0/7.0;
            frequencyArray[i] = (i+1)/(imageSize*pixelSize); // in nm^-1
            R += 1.0f;
            i++;
        }

        Plot frcPlot = new Plot("FSC curve", "Frequency (nm^-1)", "Amplitude (AU)");
        frcPlot.add( "connected circle", frequencyArray, frcCurve);
        frcPlot.setColor("red");
        frcPlot.add( "line", frequencyArray, threshold);
        frcPlot.setColor("black");
        frcPlot.addLegend("FSC curve\n1/7 threshold");
        frcPlot.setLimits(frequencyArray[0], frequencyArray[imageSize/2-1], 0.0, 1.1);
        frcPlot.show();

        Plot countPlot = new Plot("Sum count curve", "Frequency (nm^-1)", "Number of counts");
        countPlot.add( "connected circle", frequencyArray, numSumCount);
        countPlot.show();

        return 0.0d;
    }


}
