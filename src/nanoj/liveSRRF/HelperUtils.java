package nanoj.liveSRRF;


import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;

public class HelperUtils {
    private static final int MAXDEBUGVALUES = 10;

    public static float[][] floatBufferTo2D(FloatBuffer data, int dim1, int dim2) {
        float[][] outData = new float[dim1][dim2];
        for (int i = 0; i < dim1; i++)
            for (int j = 0; j < dim2; j++) {
                outData[i][j] = data.get(i * dim2 + j);
            }
        return outData;
    }

    public static FloatBuffer tofloatBuffer(float[][] float2d, int dim1,
                                            int dim2) {
        FloatBuffer outFloatBuffer = ByteBuffer.allocateDirect(4 * dim1 * dim2)
                .order(ByteOrder.nativeOrder()).asFloatBuffer();
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                outFloatBuffer.put(i * dim2 + j, float2d[i][j]);
            }
        }
        return outFloatBuffer;
    }

    public static float[] floatBufferTo1D(FloatBuffer data, int dim1) {
        float[] outData = new float[dim1];
        for (int i = 0; i < dim1; i++) {
            outData[i] = data.get(i);
        }

        return outData;
    }
    /*
     * copies data into the FloatBuffer with the data according to dim1 and dim2
     */
    public static void intofloatBuffer(float[][] data, int dim1, int dim2,
                                       FloatBuffer outFloatBuffer) {
        outFloatBuffer.rewind();
        for (int i = 0; i < dim1; i++) {
            for (int j = 0; j < dim2; j++) {
                outFloatBuffer.put(i * dim2 + j, data[i][j]);
            }
        }
        outFloatBuffer.rewind();
    }

    public static void debugWriteFloat(String path, float sumPixels,
                                       String string) {
        File file = new File(path + "testdat/" + string + "_float");

        try {
            DataOutputStream dos = new DataOutputStream(
                    new BufferedOutputStream(new FileOutputStream(file)));
            dos.writeFloat(sumPixels);
            dos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void debugWriteDouble(String path, double magMax,
                                        String string) {
        File file = new File(path + "testdat/" + string + "_double");

        try {
            DataOutputStream dos = new DataOutputStream(
                    new BufferedOutputStream(new FileOutputStream(file)));
            dos.writeDouble(magMax);
            dos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static boolean compareResults(String path, float[][] data,
                                         String string) {
        int failed = 0;
        int passed = 0;
        boolean retValue = true;
        File file = new File(path + "testdat/" + string);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(file)));
            float refvalue;
            for (int i = 0; i < data.length; i++) {
                for (int j = 0; j < data[0].length; j++) {
                    refvalue = dis.readFloat();
                    if (refvalue != data[i][j]) {
                        failed++;
                        if (failed < MAXDEBUGVALUES) {
                            System.out.println("For file " + string
                                    + " at value i = " + i + " j = " + j + ":");
                            System.out.println("Reference value id "
                                    + (i * data[0].length + j) + " is "
                                    + refvalue + " != test value of "
                                    + data[i][j]);
                        }
                        retValue = false;
                    }
                    else
                    {
                        passed++;
                    }
                }
            }

            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        if (retValue)
            System.out.println("HelperUtils.compareResults: Test " + string
                    + " passed " + passed + " of " +  data.length*data[0].length + " values and failed " + failed);
        return retValue;
    }

    public static boolean compareResults(String path, float[] data,
                                         String string) {
        int failed = 0;
        boolean retValue = true;
        File file = new File(path + "testdat/" + string);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream( new FileInputStream(file) ));
            float refvalue;
            for (int i = 0; i < data.length; i++)
            {
                refvalue = dis.readFloat();
                if (refvalue != data[i])
                {
                    if (failed++ < MAXDEBUGVALUES)
                    {
                        System.out.println("For file " + string + ":");
                        System.out.println("Reference value id " + i + " is " + refvalue + " != test value of " + data[i]);
                    }
                    retValue = false;
                }
            }

            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        if (retValue)
            System.out.println("HelperUtils.compareResults: Test " + string
                    + " passed against reference data.");
        return retValue;
    }

    public static boolean compareResults(String path, FloatBuffer data,
                                         String string) {
        int failed = 0;
        int countPassed = 0;
        int totalTests = 0;
        boolean testPassed = true;
        File file = new File(path + "testdat/" + string);
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream(file)));
            totalTests = ( int ) (file.length( )/4);
            float refvalue = 0;
            for (int i = 0; i < data.capacity(); i++)
            {
                refvalue = dis.readFloat();
                if ( refvalue != data.get(i) )
                {
                    if (failed++ < MAXDEBUGVALUES)
                    {
                        System.out.println("For file " + string + ":");
                        System.out.println("Reference value id " + i + " is " + refvalue + " != test value of " + data.get(i));
                    }
                    testPassed = false;
                } else {
                    //increment positive tests
                    countPassed++;
                }
            }

            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        //if (countPassed != totalTests)
        System.out.println("HelperUtils.compareResults: Test " + string + " passed " + countPassed + " of " + totalTests + ".");
        return testPassed;
    }

    public static boolean compareResults(String path, float data, String string)
    {
        int failed = 0;
        boolean retValue = true;
        File file = new File( path + "testdat/" + string );
        try {
            DataInputStream dis = new DataInputStream(new BufferedInputStream(
                    new FileInputStream( file ) ) );
            float refvalue = 0;

            refvalue = dis.readFloat();

            if( refvalue != data )
            {
                if( failed++ < MAXDEBUGVALUES )
                {
                    System.out.println( "For file " + string + ":");
                    System.out.println( "Reference value " + refvalue + " does not match test value of " + data );
                }
                retValue = false;
            }

            dis.close();
        } catch ( Exception e )
        {
            e.printStackTrace();
        }
        if(retValue) System.out.println("HelperUtils.compareResults: Test "
                + string + " passed against reference data.");
        return retValue;
    }

    public static void debugWriteFloatArray(String path, float[][] float2d,
                                            String string) {
        File file = new File(path + "testdat/" + string + "_float_"
                + float2d.length + "_" + float2d[0].length);
        try {
            DataOutputStream dos = new DataOutputStream(
                    new BufferedOutputStream(new FileOutputStream(file)));

            for (int i = 0; i < float2d.length; i++) {
                for (int j = 0; j < float2d[0].length; j++) {
                    dos.writeFloat(float2d[i][j]);
                }
            }

            dos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void debugWriteFloatArray(String path, float[] floatArray,
                                            String string) {
        File file = new File(path + "testdat/" + string + "_float_"
                + floatArray.length);

        try {
            DataOutputStream dos = new DataOutputStream(
                    new BufferedOutputStream(new FileOutputStream(file)));

            for (int i = 0; i < floatArray.length; i++) {
                dos.writeFloat(floatArray[i]);
            }

            dos.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void injectGoodData(FloatBuffer data, String path, String string)
    {
        File file = new File(path + "testdat/" + string);
        try {
            DataInputStream dis = new DataInputStream( new BufferedInputStream( new FileInputStream(file)) );

            for (int i = 0; i < data.capacity(); i++) {
                data.put(i, dis.readFloat() );
            }

            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void injectGoodData(float[] data, String path, String string)
    {
        File file = new File(path + "testdat/" + string);
        try {
            DataInputStream dis = new DataInputStream( new BufferedInputStream( new FileInputStream(file)) );

            for (int i = 0; i < data.length; i++) {
                data[i] = dis.readFloat();
            }

            dis.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void zeroFloatBuffer(FloatBuffer c, int capacity) {
        for (int index = 0; index < capacity; index++) {
            c.put(index, 0);
        }
    }

    public static FloatBuffer CopyFloatBuffer(FloatBuffer data) {
        FloatBuffer copy = ByteBuffer.allocateDirect(4 * data.capacity())
                .order(ByteOrder.nativeOrder()).asFloatBuffer();
        for (int i = 0; i < data.capacity(); i++) {
            copy.put(i, data.get(i));
        }

        return copy;
    }

    public static void compareResults(FloatBuffer reference, FloatBuffer test,
                                      final int diffNumber, String testName) {
        int diff = 0;
        int match = 0;
        boolean matchStatus = true;
        for (int i = 0; i < reference.capacity(); i++) {
            float a = reference.get(i);
            float b = test.get(i);
            if (a != b) {
                matchStatus = false;
                if (diff < diffNumber) {
                    System.out.println("Reference value " + i + " of " + a
                            + " != " + b + " for test " + testName);
                    diff = diff + 1;
                }
            }
            else { match++; }
        }
        if (matchStatus == true)
        {
            System.out.println("All results matched for test " + testName);
        }
        else
        {
            System.out.println("There were " + match + " total matches of " + reference.capacity() + " comparisions, for test " + testName);
        }
    }

    public static void into2DFloatBuffer(float[][] data, FloatBuffer buffer)
    {
        for(int i = 0; i < data.length; i++)
        {
            buffer.get(data[i]);
        }
    }

    // -- Convert time to string --
    public static String timeToString(double time){

        String timeString;
        int _h = (int) (time / 3600);
        int _m = (int) (((time % 86400) % 3600) / 60);
        int _s = (int) (((time % 86400) % 3600) % 60);
        if (_h > 0) timeString = _h+"h "+_m+"m "+_s+"s";
        else {
            if (_m > 0) timeString = _m+"m "+_s+"s";
            else timeString = _s+"s";
        }

        return timeString;

    }

}
