package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.IJ;

import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static java.lang.Math.min;
import static nanoj.core2.NanoJCL.*;

public class NNDistance_CL {

    // OpenCL formats
    static private CLContext context;
    static private CLProgram programNND;
    static private CLKernel kernelNND2,
    kernelSqrtNND2;

    static private CLPlatform clPlatformMaxFlop;
    static private CLDevice clDeviceMaxFlop;

    static private CLCommandQueue queue;

    private CLBuffer<FloatBuffer> clBufferXY,
            clBufferNND;


    private int nLocs,
            blockLength = 2000; // TODO: determine depending on the GPU?

    public NNDistance_CL(float[] xyDataArray){

        int nDevices = 0;
        nLocs = (xyDataArray.length)/2;
        CLPlatform[] allPlatforms;

        try{allPlatforms = CLPlatform.listCLPlatforms();}
        catch(CLException ex) {
            IJ.log("Something went wrong initializing OpenCL.");
            throw new RuntimeException("Something went wrong initializing OpenCL.");
        }

        double nFlops = 0;

        for (CLPlatform allPlatform : allPlatforms) {
            CLDevice[] allCLdeviceOnThisPlatform = allPlatform.listCLDevices();
            nDevices += allCLdeviceOnThisPlatform.length;

            for (CLDevice clDevice : allCLdeviceOnThisPlatform) {
                IJ.log("--------");
                IJ.log("Device name: " + clDevice.getName());
                IJ.log("Device type: " + clDevice.getType());
                IJ.log("Max clock: " + clDevice.getMaxClockFrequency() + " MHz");
                IJ.log("Number of compute units: " + clDevice.getMaxComputeUnits());
                IJ.log("Max work group size:: " + clDevice.getMaxWorkGroupSize());
                IJ.log("Max work item dimensions: " + clDevice.getMaxWorkItemDimensions());
                IJ.log("Max work item sizes: " + clDevice.getMaxWorkItemSizes()[0]+"/"+clDevice.getMaxWorkItemSizes()[1]+"/"+clDevice.getMaxWorkItemSizes()[2]);
                IJ.log("Global memory size: " + (clDevice.getGlobalMemSize()/1000000)+" MB");
                IJ.log("Global memory cache size: " + (clDevice.getGlobalMemCachelineSize()/1000000)+" MB");


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
        context = CLContext.create(clPlatformMaxFlop);
        CLDevice chosenDevice = context.getMaxFlopsDevice();

        clBufferXY = context.createFloatBuffer(2*nLocs, READ_ONLY);
        clBufferNND = context.createFloatBuffer(nLocs, READ_ONLY); // should always be divisable by 2

        System.out.println("Used device memory: " + (
                        clBufferXY.getCLSize() +
                        clBufferNND.getCLSize()
        )
                / 1000000d + "MB");

        String programString = getResourceAsString(liveSRRF_CL.class, "NNDistance.cl");
        programString = replaceFirst(programString, "$NLOCS$", "" + nLocs);

        programNND = context.createProgram(programString).build();
        kernelNND2 = programNND.createCLKernel("calculateNND2");
        kernelSqrtNND2 = programNND.createCLKernel("calculateSqrtNND2");

        kernelNND2.setArg(0, clBufferXY); // make sure type is the same !!
        kernelNND2.setArg(1, clBufferNND); // make sure type is the same !!
        kernelSqrtNND2.setArg(0, clBufferNND); // make sure type is the same !!

        queue = chosenDevice.createCommandQueue();

        fillBuffer(clBufferXY, xyDataArray);
        queue.putWriteBuffer(clBufferXY, false);
        queue.finish(); // Make sure everything is done
//        queue.put1DRangeKernel(kernelNND2, 0, nLocs,0);

        int workSize;
        int nBlocks = nLocs / blockLength + ((nLocs % blockLength == 0) ? 0 : 1);
        IJ.log("Number of blocks: "+nBlocks);
        System.out.println("1");
        System.out.println("Number of blocks: "+nBlocks);


        for (int nB = 0; nB < nBlocks; nB++) {
            workSize = min(blockLength, nLocs - nB * blockLength);
            queue.put1DRangeKernel(kernelNND2, nB * blockLength, workSize, 0);
            queue.finish(); // Make sure everything is done

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
            IJ.showProgress((double) nB / (double) nBlocks);
        }

        for (int nB = 0; nB < nBlocks; nB++) {
            workSize = min(blockLength, nLocs - nB * blockLength);
            queue.put1DRangeKernel(kernelSqrtNND2, nB * blockLength, workSize, 0);
            queue.finish(); // Make sure everything is done

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
            IJ.showProgress((double) nB / (double) nBlocks);
        }

    }

    public double[] readNNDbuffer() {

        queue.finish(); // Make sure everything is done
        queue.putReadBuffer(clBufferNND, true);
        FloatBuffer bufferNND = clBufferNND.getBuffer();

        double[] dataNND = new double[nLocs];
        for (int i=0; i<nLocs; i++){
            dataNND[i] = (double) bufferNND.get(i);
        }

        return dataNND;

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
