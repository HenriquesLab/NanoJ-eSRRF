package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static nanoj.core2.NanoJCL.*;

public class RemoveMacroPixelArtefacts_CL {

    static int magnification,
            widthM,
            heightM;

    public ImageStack imsResults;

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programRemoveMPartefatcs;
    static private CLKernel kernelCalculateMPmap,
            kernelCorrectMPmap;

    static private CLPlatform clPlatformMaxFlop;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer>
            clBufferImage,
            clBufferMPmap;


    // -- Check devices --
    public RemoveMacroPixelArtefacts_CL(int width, int height, int magnification) {

        this.magnification = magnification;
        this.heightM = height;
        this.widthM = width;
        this.imsResults = new ImageStack(widthM, heightM);

        CLPlatform[] allPlatforms;

        try {
            allPlatforms = CLPlatform.listCLPlatforms();
        } catch (CLException ex) {
            IJ.log("Something went wrong initializing OpenCL.");
            throw new RuntimeException("Something went wrong initializing OpenCL.");
        }

        double nFlops = 0;

        for (CLPlatform allPlatform : allPlatforms) {
            CLDevice[] allCLdeviceOnThisPlatform = allPlatform.listCLDevices();

            for (CLDevice clDevice : allCLdeviceOnThisPlatform) {
                IJ.log("--------");
                IJ.log("Device name: " + clDevice.getName());
                IJ.log("Device type: " + clDevice.getType());
                IJ.log("Max clock: " + clDevice.getMaxClockFrequency() + " MHz");
                IJ.log("Number of compute units: " + clDevice.getMaxComputeUnits());
                if (clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency() > nFlops) {
                    nFlops = clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency();
                    clPlatformMaxFlop = allPlatform;
                }
            }
        }


        context = CLContext.create(clPlatformMaxFlop);
        CLDevice chosenDevice = context.getMaxFlopsDevice();

        System.out.println("using " + chosenDevice);

        clBufferImage = context.createFloatBuffer(widthM*heightM, READ_ONLY);
        clBufferMPmap = context.createFloatBuffer(magnification*magnification, READ_WRITE);

        String programString = getResourceAsString(RemoveMacroPixelArtefacts_CL.class, "removeMPartefacts.cl");
        programString = replaceFirst(programString, "$MAGNIFICATION$", "" + magnification);
        programString = replaceFirst(programString, "$WIDTH$", "" + (widthM/magnification));
        programString = replaceFirst(programString, "$HEIGHT$", "" + (heightM/magnification));
        programString = replaceFirst(programString, "$WH$", "" + (widthM/magnification*heightM/magnification));
        programString = replaceFirst(programString, "$WM$", "" + (widthM));
        programString = replaceFirst(programString, "$HM$", "" + (heightM));
        programString = replaceFirst(programString, "$WHM$", "" + (widthM * heightM));

        programRemoveMPartefatcs = context.createProgram(programString).build();
        kernelCalculateMPmap = programRemoveMPartefatcs.createCLKernel("kernelCalculateMPmap");
        kernelCorrectMPmap = programRemoveMPartefatcs.createCLKernel("kernelCorrectMPmap");

        int argn = 0;
        kernelCalculateMPmap.setArg(argn++, clBufferImage); // make sure type is the same !!
        kernelCalculateMPmap.setArg(argn++, clBufferMPmap); // make sure type is the same !!

        argn = 0;
        kernelCorrectMPmap.setArg(argn++, clBufferImage); // make sure type is the same !!
        kernelCorrectMPmap.setArg(argn++, clBufferMPmap); // make sure type is the same !!

        queue = chosenDevice.createCommandQueue();

    }

    // --- Load the data one frame at a time ---
    public void loadData(ImageProcessor ip){

        ImageStack ims = new ImageStack(widthM, heightM);
        ims.addSlice(ip);

        fillBuffer(clBufferImage, ims);
        queue.putWriteBuffer(clBufferImage,true);

    }

    // --- Execute ---
    public void correct(){
        queue.put1DRangeKernel(kernelCalculateMPmap, 0, magnification * magnification, 0);
        queue.put1DRangeKernel(kernelCorrectMPmap, 0, widthM*heightM, 0);
        queue.finish(); // Make sure everything is done
    }

    // --- read the buffer ---
    public void readBuffer(){

        queue.putReadBuffer(clBufferImage, true);
        FloatBuffer bufferCorrectedImage = clBufferImage.getBuffer();

        float[] dataCorrectedImage;

        // Load average
        dataCorrectedImage = new float[widthM*heightM];
        for (int n = 0; n < widthM*heightM; n++) {
            dataCorrectedImage[n] = bufferCorrectedImage.get(n);
            if (Float.isNaN(dataCorrectedImage[n])) dataCorrectedImage[n] = 0; // make sure we dont get any weirdness
        }
        imsResults.addSlice(new FloatProcessor(widthM, heightM, dataCorrectedImage));


    }


    // --- Release GPU context ---
    public void release() {
        //context.release();
        while (!context.isReleased()){
            IJ.log("-------------");
            IJ.log("Releasing context...");
            context.release();
        }
    }


}
