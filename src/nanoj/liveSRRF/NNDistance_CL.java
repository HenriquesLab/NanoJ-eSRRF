package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;

import java.nio.FloatBuffer;
import java.nio.IntBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.READ_WRITE;
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.*;

public class NNDistance_CL {

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programNND;
    static private CLKernel kernelNND2,
            kernelSqrtNND2;

    static private CLPlatform clPlatformMaxFlop;
    static private CLDevice clDeviceMaxFlop,
            chosenDevice;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer> clBufferXY,
            clBufferNND;

    private CLBuffer<IntBuffer> clBufferDensity;

    private int nLocs,
            blockPerAxis,
            blockLength; // TODO: determine depending on the GPU?

    private float densityRadius, imageSize;

    private float[] xyDataArray;

    public ImageStack imsNND, imsDensity;

    public NNDistance_CL(float densityRadius, int blockPerAxis, float imageSize, int blockLength) {

        this.densityRadius = densityRadius;
        this.blockPerAxis = blockPerAxis;
        this.imageSize = imageSize;
        this.blockLength = blockLength;

        this.imsNND = new ImageStack(blockPerAxis, blockPerAxis);
        this.imsDensity = new ImageStack(blockPerAxis, blockPerAxis);


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
//                IJ.log("--------");
//                IJ.log("Device name: " + clDevice.getName());
//                IJ.log("Device type: " + clDevice.getType());
//                IJ.log("Max clock: " + clDevice.getMaxClockFrequency() + " MHz");
//                IJ.log("Number of compute units: " + clDevice.getMaxComputeUnits());
//                IJ.log("Max work group size:: " + clDevice.getMaxWorkGroupSize());
//                IJ.log("Max work item dimensions: " + clDevice.getMaxWorkItemDimensions());
//                IJ.log("Max work item sizes: " + clDevice.getMaxWorkItemSizes()[0] + "/" + clDevice.getMaxWorkItemSizes()[1] + "/" + clDevice.getMaxWorkItemSizes()[2]);
//                IJ.log("Global memory size: " + (clDevice.getGlobalMemSize() / 1000000) + " MB");
//                IJ.log("Global memory cache size: " + (clDevice.getGlobalMemCachelineSize() / 1000000) + " MB");

                if (clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency() > nFlops) {
                    nFlops = clDevice.getMaxComputeUnits() * clDevice.getMaxClockFrequency();
                    clPlatformMaxFlop = allPlatform;
                    clDeviceMaxFlop = clDevice;
                }
            }
        }

        IJ.log("--------");
        IJ.log("Maximum flops device: " + clDeviceMaxFlop.getName());

        System.out.println("Using the fastest device...");


    }

    public void calculateNND(float[] xyDataArray){

        this.xyDataArray = xyDataArray;
        context = CLContext.create(clPlatformMaxFlop);
        chosenDevice = context.getMaxFlopsDevice();

        nLocs = (xyDataArray.length) / 2;

        clBufferXY = context.createFloatBuffer(2*nLocs, READ_ONLY);
        clBufferNND = context.createFloatBuffer(nLocs, READ_WRITE);
        clBufferDensity = context.createIntBuffer(nLocs, READ_WRITE);

        System.out.println("Used device memory: " + (
                clBufferXY.getCLSize() +
                        clBufferNND.getCLSize() +
                        clBufferDensity.getCLSize()
        )
                / 1000000d + "MB");

        String programString = getResourceAsString(liveSRRF_CL.class, "NNDistance.cl");
        programString = replaceFirst(programString, "$NLOCS$", "" + nLocs);
        programString = replaceFirst(programString, "$DENSITYRADIUS2$", "" + densityRadius*densityRadius);

        programNND = context.createProgram(programString).build();
        kernelNND2 = programNND.createCLKernel("calculateNND2");
        kernelSqrtNND2 = programNND.createCLKernel("calculateSqrtNND2");

        kernelNND2.setArg(0, clBufferXY); // make sure type is the same !!
        kernelNND2.setArg(1, clBufferNND); // make sure type is the same !!
        kernelNND2.setArg(2,clBufferDensity);
        kernelSqrtNND2.setArg(0, clBufferNND); // make sure type is the same !!

        queue = chosenDevice.createCommandQueue();

        fillBuffer(clBufferXY, xyDataArray);
        queue.putWriteBuffer(clBufferXY, false);
        queue.finish(); // Make sure everything is done

        executeInBlocks(kernelNND2, blockLength, nLocs);
        executeInBlocks(kernelSqrtNND2, blockLength, nLocs);


//        int workSize;
//        int nBlocks = nLocs / blockLength + ((nLocs % blockLength == 0) ? 0 : 1);
////        IJ.log("Number of blocks: "+nBlocks);
////        System.out.println("Number of blocks: "+nBlocks);
//
//
//        for (int nB = 0; nB < nBlocks; nB++) {
//            workSize = min(blockLength, nLocs - nB * blockLength);
//            queue.put1DRangeKernel(kernelNND2, nB * blockLength, workSize, 0);
//            queue.finish(); // Make sure everything is done
//
//            if (IJ.escapePressed()) {
//                IJ.resetEscape();
//                return;
//            }
////            IJ.showProgress((double) nB / (double) nBlocks);
//        }
//
//        for (int nB = 0; nB < nBlocks; nB++) {
//            workSize = min(blockLength, nLocs - nB * blockLength);
//            queue.put1DRangeKernel(kernelSqrtNND2, nB * blockLength, workSize, 0);
//            queue.finish(); // Make sure everything is done
//
//            if (IJ.escapePressed()) {
//                IJ.resetEscape();
//                return;
//            }
//
//        }
    }

    public double[] readBuffers() {

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferNND, true);
        FloatBuffer bufferNND = clBufferNND.getBuffer();

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferDensity, true);
        IntBuffer bufferDensity = clBufferDensity.getBuffer();

        double[] dataBuffer = new double[2*nLocs];
        double[] dataNNDblocks = new double[blockPerAxis*blockPerAxis];
        double[] dataDensityBlocks = new double[blockPerAxis*blockPerAxis];
        double[] nLocsPerSpatialBlock = new double[blockPerAxis*blockPerAxis];

        int ind;

        for (int i=0; i<nLocs; i++){
            dataBuffer[i] = (double) bufferNND.get(i);
            dataBuffer[i+nLocs] = (double) bufferDensity.get(i);

            ind = (int) ((float) blockPerAxis * xyDataArray[i]/imageSize) + blockPerAxis * (int) ((float) blockPerAxis * xyDataArray[i+nLocs]/imageSize);
            nLocsPerSpatialBlock[ind] ++;
            dataNNDblocks[ind] += (double) bufferNND.get(i);
            dataDensityBlocks[ind] += (double) bufferDensity.get(i);
        }

        for (int i = 0; i<(blockPerAxis*blockPerAxis); i++){
            dataNNDblocks[i] /= nLocsPerSpatialBlock[i];
            dataDensityBlocks[i] /= nLocsPerSpatialBlock[i];
        }

        imsNND.addSlice(new FloatProcessor(blockPerAxis, blockPerAxis, dataNNDblocks));
        imsDensity.addSlice(new FloatProcessor(blockPerAxis, blockPerAxis, dataDensityBlocks));

        return dataBuffer;
    }


    // --- Release GPU context ---
    public void release() {
        //context.release();
        while (!context.isReleased()){
//            IJ.log("-------------");
//            IJ.log("Releasing context...");
            context.release();
        }
    }


    private void executeInBlocks(CLKernel kernel, int blockSize, int totalNexecutions){

        int workSize;
        int nBlocks = totalNexecutions / blockSize + ((totalNexecutions % blockSize == 0) ? 0 : 1);

        for (int nB = 0; nB < nBlocks; nB++) {
            workSize = min(blockSize, totalNexecutions - nB * blockSize);
            queue.put1DRangeKernel(kernel, nB * blockSize, workSize, 0);
            queue.finish(); // Make sure everything is done

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }

    }


}
