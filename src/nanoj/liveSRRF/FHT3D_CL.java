package nanoj.liveSRRF;

import java.nio.FloatBuffer;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;

import com.jogamp.opencl.CLBuffer;
import com.jogamp.opencl.CLCommandQueue;
import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLPlatform;
import com.jogamp.opencl.CLDevice.Type;
import com.jogamp.opencl.CLKernel;
import com.jogamp.opencl.CLProgram;
import com.jogamp.opencl.CLMemory.Mem;

import ij.IJ;
import ij.ImageStack;
import ij.process.FloatProcessor;


// Original 3D FHT code from LOCI: OpenCL implementation of 3D FHT for the purpose of deconvolution
// Code extracted from https://github.com/uw-loci/opencl-decon/blob/master/src/demos/FHT3DExample.java
// Re-implemented here by Romain Laine 2019-12-16, r.laine@ucl.ac.uk

public class FHT3D_CL {

    private static boolean DEBUG = false;

    private CLContext context;
    private CLBuffer<FloatBuffer> clFloatBufferDATA;
    private CLBuffer<FloatBuffer> clFloatBufferS;
    private CLBuffer<FloatBuffer> clFloatBufferC;
    private CLBuffer<FloatBuffer> clFloatBufferSW;
    private CLBuffer<FloatBuffer> clFloatBufferCW;
    private CLBuffer<FloatBuffer> clFloatBufferSH;
    private CLBuffer<FloatBuffer> clFloatBufferCH;
    private CLProgram programFHT3D;
    private CLKernel kernelsw;
    private CLKernel kernelsh;
    private CLKernel kernels;
    private CLKernel kernelf;
    private CLCommandQueue queue;
    private int globalWorkSize = 1;
    private int localWorkSize = 1;
    private int maxComputeUnits = 0;
    private FloatBuffer dataA;
    private int w = 0;
    private int h = 0;
    private int d = 0;
    private int inverse = 0;

    public static ImageStack imsFHT;

    /*
     * Constructor used for OpenCL persistence, sets up GPU context and allocates arrays based on size of recent calls
     *
     * it takes in the properties of the inputs but does not calculate the results.
     *
     */

    protected void release() {
        dataA = null;
        queue.release();
    }

    //protected void finalize() throws Throwable
    //{
    //	dataA = null;
    //	queue.release();
    //}
    public static boolean powerOf2Size(int w) {
        {
            int i = 2;
            while (i < w) {
                i = i * 2;
            }
            return i == w;
        }
    }

    public FHT3D_CL(int w, int h, int d) {
        if (DEBUG) {
            System.out.println("Asserting that dimensions meet design criteria");
        }
        assert (w == h);
        assert (h == d);
        assert (powerOf2Size(d));
        assert (d <= 1024);

        this.w = w;
        this.h = h;
        this.d = d;

        imsFHT = new ImageStack(w, h);

        try {
            if (DEBUG) {
                System.out.println("Local work size dimensions are max array size of");
            }
            CLPlatform[] platforms = CLPlatform.listCLPlatforms();
            if (DEBUG) {
                for (CLPlatform clPlatform : platforms) {
                    IJ.log("Discovered " + clPlatform.getName());
                }
            }
            context = CLContext.create(Type.GPU);

//            InputStream programStream = FHT3D_CL.class.getResourceAsStream("fht.cl"); // TODO: somehow the CL file cannot be found when in the resources folder
//            IJ.log("PROGRAM: "+programStream);
            this.programFHT3D = context.createProgram(FHT3D_CL.class.getResourceAsStream("fht.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (DEBUG) {
            System.out.println("assign the kernel");
        }
        this.kernelsw = programFHT3D.createCLKernel("fhtsw");
        this.kernelsh = programFHT3D.createCLKernel("fhtsh");
        this.kernels = programFHT3D.createCLKernel("fhts");
        this.kernelf = programFHT3D.createCLKernel("fhtf");

        if (DEBUG) {
            System.out.println("get a command queue");
        }
        this.queue = context.getMaxFlopsDevice().createCommandQueue(CLCommandQueue.Mode.PROFILING_MODE);

        if (DEBUG) {
            System.out.println("find the max compute unites");
        }
        maxComputeUnits = queue.getDevice().getMaxComputeUnits();

        this.dataA = ByteBuffer.allocateDirect(d * w * h * 4).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer sw = ByteBuffer.allocateDirect(w).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer cw = ByteBuffer.allocateDirect(w).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer sh = ByteBuffer.allocateDirect(h).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer ch = ByteBuffer.allocateDirect(h).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer s = ByteBuffer.allocateDirect(d).order(ByteOrder.nativeOrder()).asFloatBuffer();
        FloatBuffer c = ByteBuffer.allocateDirect(d).order(ByteOrder.nativeOrder()).asFloatBuffer();

        if (DEBUG) {
            System.out.println("Generating sin cos tables");
        }
        makeSinCosTables(d, s, c);
        makeSinCosTables(w, sw, cw);
        makeSinCosTables(h, sh, ch);

        if (DEBUG) {
            System.out.println("Create a CLMemory Structure for  c, w, ...");
        }
        clFloatBufferC = context.createBuffer(c, Mem.READ_ONLY);
        clFloatBufferS = context.createBuffer(s, Mem.READ_ONLY);
        clFloatBufferCW = context.createBuffer(cw, Mem.READ_ONLY);
        clFloatBufferSW = context.createBuffer(sw, Mem.READ_ONLY);
        clFloatBufferCH = context.createBuffer(ch, Mem.READ_ONLY);
        clFloatBufferSH = context.createBuffer(sh, Mem.READ_ONLY);
        clFloatBufferDATA = context.createFloatBuffer(d * w * h, Mem.READ_WRITE);

        if (DEBUG) {
            System.out.println("Making kernel assignments");
        }
        kernels.setArg(0, clFloatBufferDATA);
        kernelsw.setArg(0, clFloatBufferDATA);
        kernelsh.setArg(0, clFloatBufferDATA);
        kernelf.setArg(0, clFloatBufferDATA);
        kernels.setArg(1, clFloatBufferS);
        kernels.setArg(2, clFloatBufferC);
        kernelsw.setArg(1, clFloatBufferSW);
        kernelsw.setArg(2, clFloatBufferCW);
        kernelsh.setArg(1, clFloatBufferSH);
        kernelsh.setArg(2, clFloatBufferCH);
        kernelsw.setArg(3, w);
        kernelsw.setArg(4, h);
        kernelsw.setArg(5, d);
        kernels.setArg(3, w);
        kernels.setArg(4, h);
        kernels.setArg(5, d);
        kernelsh.setArg(3, w);
        kernelsh.setArg(4, h);
        kernelsh.setArg(5, d);
        kernelf.setArg(1, w);
        kernelf.setArg(2, h);
        kernelf.setArg(3, d);
        kernels.setArg(6, maxComputeUnits);
        kernelsw.setArg(6, maxComputeUnits);
        kernelsh.setArg(6, maxComputeUnits);
        kernelf.setArg(4, maxComputeUnits);

        kernelsh.setNullArg(7, h * 4);
        kernelsh.setNullArg(8, h);
        kernelsh.setNullArg(9, h);
        kernelsw.setNullArg(7, w * 4);
        kernelsw.setNullArg(8, w);
        kernelsw.setNullArg(9, w);
        kernels.setNullArg(7, d * 4);
        kernels.setNullArg(8, d);
        kernels.setNullArg(9, d);

        if (DEBUG) {
            System.out.println("Copying data to the device");
        }
        queue.putWriteBuffer(clFloatBufferC, false);
        queue.putWriteBuffer(clFloatBufferS, false);
        queue.putWriteBuffer(clFloatBufferCW, false);
        queue.putWriteBuffer(clFloatBufferSW, false);
        queue.putWriteBuffer(clFloatBufferCH, false);
        queue.putWriteBuffer(clFloatBufferSH, false);
    }


    // --- Run, Forest, run ---
    public synchronized long run(ImageStack ims, boolean inverse) {

        if (DEBUG) {
            System.out.println("add the input data to the objects buffer ");
        }
        for (int z = 0; z < d; z++) {
            for (int y = 0; y < h; y++) {
                for (int x = 0; x < w; x++) {
                    dataA.put(z * w * h + y * w + x, (float) ims.getVoxel(x,y,z));
                }
            }
        }

        long executionTime = run(dataA, inverse);

        return executionTime;

    }


    public synchronized long run(float[][] data, boolean inverse) {

        if (DEBUG) {
            System.out.println("add the input data to the objects buffer ");
        }
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < w * h; j++) {
                dataA.put(i * w * h + j, data[i][j]);
            }
        }

        long executionTime = run(dataA, inverse);


        if (DEBUG) {
            System.out.println("Copy data to the float[][]...");
        }
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < w * h; j++) {
                data[i][j] = dataA.get(i * w * h + j);
            }
        }

        return executionTime;

    }

    public synchronized long run(FloatBuffer data, boolean inverse) {
        clFloatBufferDATA.use(data);

        if (DEBUG) {
            System.out.println("translate for inverse");
        }
        this.inverse = (inverse == true) ? 1 : 0;

        if (DEBUG) {
            System.out.println("Setting CLBuffer to input data buffer...");
        }
        queue.putWriteBuffer(clFloatBufferDATA, false);

        long executionTime = run();

        if (DEBUG) {
            System.out.println("Reading back data from the device...");
        }
        queue.putReadBuffer(clFloatBufferDATA, true);
        FloatBuffer bufferFHT = clFloatBufferDATA.getBuffer();

        for (int f = 0; f < d; f++) {
            float[] dataFHT = new float[w * h];
            for (int i = 0; i < w*h; i++) {
                dataFHT[i] = bufferFHT.get(i + f*w*h);
//                if (Float.isNaN(dataFHT[i])) dataFHT[i] = 0; // make sure we dont get any weirdness
            }
            imsFHT.addSlice(new FloatProcessor(w, h, dataFHT));
        }

        return executionTime;
    }

    private synchronized long run() {

        //this changes every time
        kernelf.setArg(5, inverse);

        if (DEBUG) {
            System.out.println("Setting global max and mins...");
        }
        globalWorkSize = w * maxComputeUnits;
        localWorkSize = w;

        if (DEBUG) {
            System.out.println("Writing the input data set to the device...");
        }
        //put the data on the card
//        queue.putWriteBuffer(clFloatBufferDATA, false); // duplicate line with above run method.

        long executionTime = System.currentTimeMillis();

        if (DEBUG) {
            System.out.println("Enqueing the kernel w...");
        }
        queue.put1DRangeKernel(kernelsw, 0, globalWorkSize, localWorkSize); // TODO: split worksizes to fit on GPU?

        //if( DEBUG ) { System.out.println("Setting global max and mins...");}
        //globalWorkSize = h*maxComputeUnits;
        //localWorkSize = h;

        if (DEBUG) {
            System.out.println("Enqueing the kernel h...");
        }
        queue.put1DRangeKernel(kernelsh, 0, globalWorkSize, localWorkSize);

        //if( DEBUG ) { System.out.println("Setting global max and mins...");}
        //globalWorkSize = d*maxComputeUnits;
        //localWorkSize = d;

        if (DEBUG) {
            System.out.println("Enqueing the kernel s...");
        }
        queue.put1DRangeKernel(kernels, 0, globalWorkSize, localWorkSize);

        if (DEBUG) {
            System.out.println("Enqueing the kernel f...");
        }
        queue.put1DRangeKernel(kernelf, 0, globalWorkSize, localWorkSize);

        if (DEBUG) {
            System.out.println("Waiting for CLFinish to return...");
        }
        queue.finish();

        executionTime = System.currentTimeMillis() - executionTime;

        return executionTime;
    }

//    public double getExecutionTimeInNS(CLEvent clevent) { // not used in this code at teh moment
//        long start;
//        long end;
//
//        end = clevent.getProfilingInfo(ProfilingCommand.END);
//        start = clevent.getProfilingInfo(ProfilingCommand.START);
//
//        // convert nanoseconds to seconds on return
//        return (double) 1.0e-9 * (end - start);
//    }


    public static void makeSinCosTables(int maxN, FloatBuffer s, FloatBuffer c) {
        int n = maxN / 4;
        double theta = 0.0;
        double dTheta = 2.0 * Math.PI / maxN;
        for (int i = 0; i < n; i++) {
            c.put(i, (float) Math.cos(theta));
            s.put(i, (float) Math.sin(theta));
            theta += dTheta;
        }
    }
}
