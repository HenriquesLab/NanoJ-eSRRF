package nanoj.liveSRRF;

public class LinearRegressions {

    public static double linearRegressionLeastSquareNoOffset(double[] xArray, double[] yArray){
        // see info on http://mathworld.wolfram.com/LeastSquaresFitting.html
        // this linear regression assumes no offset

        double XX = 0;
        double XY = 0;

        for (int i = 0; i < xArray.length; i++) {
            XX += xArray[i]*xArray[i];
            XY += xArray[i]*yArray[i];
        }
        return XY/XX; // this returns the slope estimated from the least square method
    }

    public static double linearRegressionLeastSquareNoOffset(double[] yArray){
        // see info on http://mathworld.wolfram.com/LeastSquaresFitting.html
        // this linear regression assumes no offset (crosses the 0,0 origin of xes)

        // Find the index with the closest value to zero!
        double dZero = Double.MAX_VALUE;
        int idZero = 0;
        for (int i = 0; i < yArray.length; i++) {
            if (Math.abs(yArray[i]) < dZero){
                dZero = Math.abs(yArray[i]);
                idZero = i;
            }
        }

        // Build an xArray that crosses 0 in the right place for the fitting
        double[] xArray = new double[yArray.length];
        for (int i = 0; i < yArray.length; i++) {
            xArray[i] = i - idZero;
//            IJ.log("xArray: "+xArray[i]);
        }

        double XX = 0;
        double XY = 0;

        for (int i = 0; i < xArray.length; i++) {
            XX += xArray[i]*xArray[i];
            XY += xArray[i]*yArray[i];
        }
        return XY/XX; // this returns the slope estimated from the least square method
    }


    public static float linearRegressionLeastSquareNoOffset(float[] xArray, float[] yArray){
        // see info on http://mathworld.wolfram.com/LeastSquaresFitting.html
        // this linear regression assumes no offset

        float XX = 0;
        float XY = 0;

        for (int i = 0; i < xArray.length; i++) {
            XX += xArray[i]*xArray[i];
            XY += xArray[i]*yArray[i];
        }
        return XY/XX; // this returns the slope estimated from the least square method
    }

    public static float[] linearRegressionLeastSquare(float[] xArray, float[] yArray){
        // see info on http://mathworld.wolfram.com/LeastSquaresFitting.html

        int n = xArray.length;
        float x_bar = 0;
        float y_bar = 0;
        float x2 = 0;
        float xy = 0;

        for (int i=0; i<n; i++){
            x_bar += xArray[i]/(float) n;
            y_bar += yArray[i]/(float) n;
            x2 += xArray[i]*xArray[i];
            xy += xArray[i]*yArray[i];
        }
        // returns {offset, slope}
        return new float[]{(y_bar * x2 - x_bar * xy) / (x2 - n * x_bar * x_bar), (xy - n * x_bar * y_bar) / (x2 - n * x_bar * x_bar)};
    }

    public static double[] linearRegressionLeastSquare(double[] xArray, double[] yArray){
        // see info on http://mathworld.wolfram.com/LeastSquaresFitting.html

        int n = xArray.length;
        double x_bar = 0;
        double y_bar = 0;
        double x2 = 0;
        double xy = 0;

        for (int i=0; i<n; i++){
            x_bar += xArray[i]/(float) n;
            y_bar += yArray[i]/(float) n;
            x2 += xArray[i]*xArray[i];
            xy += xArray[i]*yArray[i];
        }
        // returns {offset, slope}
        return new double[]{(y_bar * x2 - x_bar * xy) / (x2 - n * x_bar * x_bar), (xy - n * x_bar * y_bar) / (x2 - n * x_bar * x_bar)};
    }



}
