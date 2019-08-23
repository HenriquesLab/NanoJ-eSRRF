package nanoj.liveSRRF;

import ij.ImageStack;
import ij.process.FloatProcessor;

// This class computes 2D stack projections on the CPU by accessing directly the float values of an imageStack.
// Let's populate the capabilities of this accordingly.
public class StackProjections {

    // Initialisation method, currently does nothing
//    public void StackProjections(){}

    // Average projections
    public static FloatProcessor calculateAverage(ImageStack ims){

        float[] avgArray = new float[ims.getWidth()*ims.getHeight()];
        int nSlices = ims.getSize();
        float[] pixels = new float[nSlices];

        for (int y=0; y<ims.getHeight(); y++){
            for (int x=0; x<ims.getWidth(); x++){
                ims.getVoxels(x, y, 0, 1,1, nSlices, pixels); // pulls the array of data from a single x,y position over the stack
                for (int s=0; s<nSlices; s++){
                    avgArray[x + y*ims.getHeight()] += pixels[s]/nSlices;
                }
            }
        }
        return new FloatProcessor(ims.getWidth(), ims.getHeight(), avgArray);
    }


    // Standard deviation projections // TODO: untested as yet
    public static FloatProcessor calculateStandardDeviation(ImageStack ims){

        double avg, avgSquare;
        float[] stdArray = new float[ims.getWidth()*ims.getHeight()];

        int nSlices = ims.getSize();
        float[] pixels = new float[nSlices];

        for (int x=0; x<ims.getWidth(); x++){
            for (int y=0; y<ims.getHeight(); y++){
                avg = 0;
                avgSquare = 0;
                ims.getVoxels(x, y, 0, 1,1, nSlices, pixels); // pulls the array of data from a single x,y position over the stack
                for (int s=0; s<nSlices; s++){
                    avg += pixels[s]/nSlices;
                    avgSquare += pixels[s]*pixels[s]/nSlices;
                }
                stdArray[x + y*ims.getHeight()] = (float) Math.sqrt(avgSquare - avg*avg);
            }
        }
        return new FloatProcessor(ims.getWidth(), ims.getHeight(), stdArray);
    }

    // Maximum intensity projection
    public static FloatProcessor calculateMIP(ImageStack ims){

        float mip;
        float[] mipArray = new float[ims.getWidth()*ims.getHeight()];

        int nSlices = ims.getSize();
        float[] pixels = new float[nSlices];

        for (int x=0; x<ims.getWidth(); x++){
            for (int y=0; y<ims.getHeight(); y++){
                mip = 0;
                ims.getVoxels(x, y, 0, 1,1, nSlices, pixels); // pulls the array of data from a single x,y position over the stack
                for (int s=0; s<nSlices; s++){
                    if (pixels[s]>mip) mip = pixels[s];
                }
                mipArray[x + y*ims.getHeight()] = mip;
            }
        }
        return new FloatProcessor(ims.getWidth(), ims.getHeight(), mipArray);
    }

}
