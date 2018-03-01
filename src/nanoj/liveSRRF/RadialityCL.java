package nanoj.liveSRRF;

import com.jogamp.opencl.*;
import ij.ImagePlus;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.io.IOException;
import java.nio.FloatBuffer;

import static com.jogamp.opencl.CLMemory.Mem.READ_ONLY;
import static com.jogamp.opencl.CLMemory.Mem.WRITE_ONLY;
import static java.lang.Math.*;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.System.nanoTime;

public class RadialityCL {

    public boolean DEBUG = false;

    static CLContext context;
    static CLProgram program;
    static CLKernel kernelCalculateRadiality;
    static CLCommandQueue queue;

    private float[] xRingCoordinates, yRingCoordinates;

    private final int width, height, widthM, heightM, magnification;
    private final double ringRadius, radialitySensitivity;

    public int nVectors = 12;

    private CLBuffer<FloatBuffer>
            clBufferPx,
            clBufferGx, clBufferGy,
            clBufferXRingCoordinates, clBufferYRingCoordinates,
            clBufferRad, clBufferIntpolated;

    public RadialityCL(int width, int height, int magnification, double ringRadius, double radialitySensitivity) {

        this.width = width;
        this.height = height;
        this.widthM = width * magnification;
        this.heightM = height * magnification;
        this.magnification = magnification;
        this.ringRadius = ringRadius;
        this.radialitySensitivity = radialitySensitivity;

        //Create a context from GPU
        //context = CLContext.create( CLDevice.Type.GPU );
        context = CLContext.create();
        System.out.println("created " + context);

        // Select fastest device
        CLDevice device = context.getMaxFlopsDevice();
        System.out.println("using " + device);

        // create command queue on device.
        queue = device.createCommandQueue();

        // create the program
        try {
            program = context.createProgram(RadialityCL.class.getResourceAsStream("/Radiality.cl")).build();
        } catch (IOException e) {
            e.printStackTrace();
        }
        kernelCalculateRadiality = program.createCLKernel("calculateRGC");

        clBufferXRingCoordinates = context.createFloatBuffer(nVectors, READ_ONLY);
        clBufferYRingCoordinates = context.createFloatBuffer(nVectors, READ_ONLY);
        clBufferPx = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferGx = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferGy = context.createFloatBuffer(width * height, READ_ONLY);
        clBufferRad = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);
        clBufferIntpolated = context.createFloatBuffer(widthM * heightM, WRITE_ONLY);

        System.out.println("used device memory: " + (
                clBufferXRingCoordinates.getCLSize() +
                        clBufferYRingCoordinates.getCLSize() +
                        clBufferPx.getCLSize() +
                        clBufferGx.getCLSize() +
                        clBufferGy.getCLSize() +
                        clBufferRad.getCLSize() +
                        clBufferIntpolated.getCLSize())
                / 1000000d + "MB");

        // pre-prepare the Ring-Coordinates
        buildRing(ringRadius*magnification);
    }

    public synchronized FloatProcessor[] calculateRadiality(ImageProcessor ip, float shiftX, float shiftY) {
        assert (ip.getWidth() == width && ip.getHeight() == height);

        FloatProcessor fpPx = ip.convertToFloatProcessor();
        FloatProcessor fpGx = new FloatProcessor(width, height);
        FloatProcessor fpGy = new FloatProcessor(width, height);
        FloatProcessor fpRadiality = new FloatProcessor(widthM, heightM);
        FloatProcessor fpInterpolated = new FloatProcessor(widthM, heightM);

        // calculate gradient
        int width1 = width-1;
        int height1 = height-1;
        for (int y=0; y<height; y++) {
            for (int x=0; x<width; x++) {
                int x0 = max(x-1, 0);
                int x1 = min(x+1, width1);
                int y0 = max(y-1, 0);
                int y1 = min(y+1, height1);
                fpGx.setf(x,y, -fpPx.getf(x0, y)+fpPx.getf(x1,y));
                fpGy.setf(x,y, -fpPx.getf(x, y0)+fpPx.getf(x,y1));
            }
        }

        if (DEBUG) {
            new ImagePlus("Gx", fpGx).show();
            new ImagePlus("Gy", fpGy).show();
        }

        // prepare CL buffers
        fillBuffer(clBufferXRingCoordinates, xRingCoordinates);
        fillBuffer(clBufferYRingCoordinates, yRingCoordinates);
        fillBuffer(clBufferPx, fpPx);
        fillBuffer(clBufferGx, fpGx);
        fillBuffer(clBufferGy, fpGy);

        // make kernelCalculateRGC assignment
        int n = 0;
        kernelCalculateRadiality.setArg( n++, clBufferXRingCoordinates ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferYRingCoordinates ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferPx ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferGx ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferGy ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferRad ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, clBufferIntpolated); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, magnification ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, (float) (ringRadius*magnification) ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, (float) (radialitySensitivity) ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, shiftX ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, shiftY ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, nVectors ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, width ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, height ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, widthM ); // make sure type is the same !!
        kernelCalculateRadiality.setArg( n++, heightM ); // make sure type is the same !!

        // asynchronous write of data to GPU device,
        // followed by blocking read to get the computed results back.
        long time = nanoTime();
        queue.putWriteBuffer( clBufferXRingCoordinates, false );
        queue.putWriteBuffer( clBufferYRingCoordinates, false );
        queue.putWriteBuffer( clBufferPx, false );
        queue.putWriteBuffer( clBufferGx, false );
        queue.putWriteBuffer( clBufferGy, false );
        queue.put2DRangeKernel(kernelCalculateRadiality, 0, 0, widthM, heightM, 0, 0);
        queue.finish();
        queue.putReadBuffer(clBufferRad, true);
        queue.putReadBuffer(clBufferIntpolated, true);
        time = nanoTime() - time;

        //clBufferRad.getBuffer().get((float[]) fpRadiality.getPixels());
        grabBuffer(clBufferRad, fpRadiality);
        grabBuffer(clBufferIntpolated, fpInterpolated);
        return new FloatProcessor[] {fpRadiality, fpInterpolated};
    }

    public void release() {
        context.release();
    }

    private void buildRing(double spatialRadius) {
        xRingCoordinates = new float[nVectors];
        yRingCoordinates = new float[nVectors];
        double angleStep = (PI * 2d) / nVectors;
        for(int angleIter = 0; angleIter < nVectors; angleIter++){
            xRingCoordinates[angleIter] = (float) (spatialRadius * cos(angleStep * angleIter));
            yRingCoordinates[angleIter] = (float) (spatialRadius * sin(angleStep * angleIter));
        }
    }

    public static void fillBuffer(CLBuffer<FloatBuffer> clBuffer, float[] data) {
        FloatBuffer buffer = clBuffer.getBuffer();
        for(int n=0; n<data.length; n++) buffer.put(n, data[n]);
    }

    public static void fillBuffer(CLBuffer<FloatBuffer> clBuffer, ImageProcessor ip) {
        FloatBuffer buffer = clBuffer.getBuffer();
        for(int n=0; n<ip.getPixelCount(); n++) buffer.put(n, ip.getf(n));
    }

    public static void grabBuffer(CLBuffer<FloatBuffer> clBuffer, FloatProcessor fp) {
        FloatBuffer buffer = clBuffer.getBuffer();
        for(int n=0; n<fp.getPixelCount(); n++) fp.setf(n, buffer.get(n));
    }
}
