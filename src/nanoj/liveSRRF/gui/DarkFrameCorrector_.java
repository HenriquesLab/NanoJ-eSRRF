package nanoj.liveSRRF.gui;

import ij.*;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import ij.process.FloatProcessor;
import nanoj.core.java.projections.Projections2D;
import nanoj.liveSRRF.RadialityCL;


/**
 * Created by sculley on 26/02/2018.
 */
public class DarkFrameCorrector_ implements PlugIn {

    Prefs prefs = new Prefs();
    String prefsHeader = this.getClass().getName();

    String titleDark, titleData;
    int magnification;
    double ringRadius, sensitivity;
    boolean doAverage, doMaximum, doStd;

    ImagePlus impDark, impData;
    Projections2D projector = new Projections2D();

    @Override
    public void run(String arg) {

        String[] imageTitles = WindowManager.getImageTitles();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Quick and dirty sCMOS correction");

        gd.addChoice("Dark frames image", imageTitles, getDefaultImageCheck("titleDark", imageTitles));

        gd.addChoice("Data frames image", imageTitles, getDefaultImageCheck("titleData", imageTitles));

        gd.addMessage("Radiality settings");
        gd.addNumericField("Magnification", getPrefs("magnification", 5),0);
        gd.addNumericField("Ring radius", getPrefs("ringRadius", 1.0f), 1);
        gd.addNumericField("Radiality sensitivity", getPrefs("sensitivity", 10), 0);

        gd.addMessage("Radiality projection options");
        gd.addCheckbox("Average", getPrefs("averageRadiality", true));
        gd.addCheckbox("Maximum Intensity", getPrefs("maximumRadiality", true));
        gd.addCheckbox("Standard Deviation", getPrefs("stdRadiality", true));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        titleDark = gd.getNextChoice();
        setPrefs("titleDark", titleDark);

        titleData = gd.getNextChoice();
        setPrefs("titleData", titleData);

        magnification = (int) gd.getNextNumber();
        setPrefs("magnification", (float) magnification);
        ringRadius = gd.getNextNumber();
        setPrefs("ringRadius", (float) ringRadius);
        sensitivity = gd.getNextNumber();
        setPrefs("sensitivity", (float) sensitivity);

        doAverage = gd.getNextBoolean();
        setPrefs("doAverage", doAverage);
        doMaximum = gd.getNextBoolean();
        setPrefs("doMaximum", doMaximum);
        doStd = gd.getNextBoolean();
        setPrefs("doStd", doStd);

        prefs.savePreferences();

        // set-up
        impDark = WindowManager.getImage(titleDark);
        impData = WindowManager.getImage(titleData);

        ImageStack imsDark = impDark.getImageStack().convertToFloat();
        ImageStack imsData = impData.getImageStack().convertToFloat();

        int wDark = imsDark.getWidth();
        int hDark = imsDark.getHeight();
        int nSlicesDark = imsDark.getSize();

        int wData = imsData.getWidth();
        int hData = imsData.getHeight();
        int nSlicesData = imsData.getSize();

        // check image dimensions
        assert (wDark==wData && hDark==hData);

        // array with projections options
        boolean[] projections = new boolean[]{doAverage, doMaximum, doStd};

        int nProjections = 0;
        for(int i=0; i<projections.length; i++){
            if(projections[i]) nProjections++;
        }

        // perform radiality transform on dark frames
        ImageStack imsDarkRadiality = new ImageStack(wDark*magnification, hDark*magnification);
        runRadiality(imsDark, imsDarkRadiality);

        IJ.log("Dark radiality done!");

        // perform dark raw data subtraction
        FloatProcessor fpDarkAverage = averageProjection(imsDark);
        ImageStack imsSubtracted = imsData.duplicate();
        subtract(imsSubtracted, fpDarkAverage);

        // average project dark radiality and interpolated intensity
        FloatProcessor fpDarkRadiality = averageProjection(imsDarkRadiality);

        // perform radiality on data frames
        ImageStack imsDataRadiality = new ImageStack(wData*magnification, hData*magnification);
        ImageStack imsDataInterpolated = new ImageStack(wData*magnification, hData*magnification);
        runRadiality(imsSubtracted, imsDataRadiality, imsDataInterpolated);
        IJ.log("Data radiality done!");

        // generate reference uncorrected images
        ImageStack imsProjectionsUncorrected = new ImageStack(wData*magnification, hData*magnification, nProjections);
        doZProjections(imsProjectionsUncorrected, imsDataRadiality);
        ImagePlus impUncorrected = new ImagePlus("Uncorrected SRRF projections", imsProjectionsUncorrected);

        // subtract average dark frame radiality from each frame
        subtract(imsDataRadiality, fpDarkRadiality);

        // perform intensity weighting
        ImageStack imsWeightedRadiality = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);

        // project data
        ImageStack imsProjections = new ImageStack(wData*magnification, hData*magnification, nProjections);
        doZProjections(imsProjections, imsWeightedRadiality);
        ImagePlus impCorrected = new ImagePlus("Corrected images", imsProjections);

        impCorrected.show();
        impUncorrected.show();

    }


    // helper function for title handling

    public String getDefaultImageCheck(String key, String[] titles){
        String targetImage = getPrefs(key, titles[0]);
        for(String s:titles){
            if(s==targetImage) return targetImage;
        }
        return titles[0];
    }

    //helper functions for prefs handling

    public float getPrefs(String key, float defaultValue) {
        return (float) prefs.get(prefsHeader+"."+key, defaultValue);
    }

    public boolean getPrefs(String key, boolean defaultValue) {
        return prefs.get(prefsHeader+"."+key, defaultValue);
    }

    public String getPrefs(String key, String defaultValue) { return prefs.get(prefsHeader+"."+key, defaultValue); }

    public void setPrefs(String key, float defaultValue) {
        prefs.set(prefsHeader+"."+key, defaultValue);
    }

    public void setPrefs(String key, boolean defaultValue) {
        prefs.set(prefsHeader+"."+key, defaultValue);
    }

    public void setPrefs(String key, String defaultValue) {prefs.set(prefsHeader+"."+key, defaultValue); }

    // helper functions for image maths

    public ImageStack calculateIntensityWeighting(ImageStack imsRadiality, ImageStack imsInterpolated){
        int nSlices = imsRadiality.getSize();
        int w = imsRadiality.getWidth();
        int h = imsRadiality.getHeight();

        ImageStack imsIW = new ImageStack(w, h, nSlices);

        for(int n=1; n<=nSlices; n++){
            float[] radialityPixels = (float[]) imsRadiality.getProcessor(n).convertToFloatProcessor().getPixels();
            float[] interpolatedPixels = (float[]) imsInterpolated.getProcessor(n).convertToFloatProcessor().getPixels();
            float[] IWPixels = new float[w*h];
            for(int i=0; i<w*h; i++){
                IWPixels[i] = radialityPixels[i]*interpolatedPixels[i];
            }
            imsIW.setProcessor(new FloatProcessor(w, h, IWPixels), n);
        }
        return imsIW;
    }

    public void subtract(ImageStack ims, FloatProcessor fp){
        int nSlices = ims.getSize();
        int w = ims.getWidth();
        int h = ims.getHeight();

        float[] fpPixels = (float[]) fp.getPixels();

        for(int n=1; n<=nSlices; n++){
            float[] imsPixels = (float[]) ims.getProcessor(n).convertToFloatProcessor().getPixels();

            for(int i=0; i<w*h; i++){
                imsPixels[i] -= fpPixels[i];
            }

            ims.setProcessor(new FloatProcessor(w,h,imsPixels), n);
        }
    }

    // helper functions for z-projecting

    public FloatProcessor averageProjection(ImageStack ims){
        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        float[] average = new float[w*h];

        for(int n=1; n<=nSlices; n++){
            FloatProcessor fp = ims.getProcessor(n).convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();
            for(int i=0; i<w*h; i++){
                average[i]+= pixels[i]/nSlices;
            }
        }
        return new FloatProcessor(w, h, average);
    }

    public FloatProcessor maxProjection(ImageStack ims){
        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        float[] maxIntensity = new float[w*h];

        for(int n=1; n<=nSlices; n++){
            FloatProcessor fp = ims.getProcessor(n).convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();
            for(int i=0; i<w*h; i++){
                maxIntensity[i] = Math.max(maxIntensity[i], pixels[i]);
            }
        }
        return new FloatProcessor(w, h, maxIntensity);

    }

    public FloatProcessor stdProjection(ImageStack ims){

        FloatProcessor fpAverage = averageProjection(ims);
        float[] average = (float[]) fpAverage.getPixels();

        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        float[] std = new float[w*h];

        for(int n=1; n<=nSlices; n++){
            FloatProcessor fp = ims.getProcessor(n).convertToFloatProcessor();
            float[] pixels = (float[]) fp.getPixels();
            for(int i=0; i<w*h; i++){
                float diff = pixels[i] - average[i];
                std[i]+= (diff*diff)/nSlices;
            }
        }

        for(int i=0; i<w*h; i++) std[i] = (float) Math.sqrt(std[i]);

        return new FloatProcessor(w, h, std);

    }

    public void doZProjections(ImageStack imsProjections, ImageStack imsRaw){

        FloatProcessor fp;
        int n=1;

        if (doAverage) {
            fp = averageProjection(imsRaw);
            imsProjections.setSliceLabel("Average", n);
            imsProjections.setProcessor(fp, n++);
        }
        if (doMaximum) {
            fp = maxProjection(imsRaw);
            imsProjections.setSliceLabel("Maximum", n);
            imsProjections.setProcessor(fp, n++);
        }
        if (doStd) {
            fp = stdProjection(imsRaw);
            imsProjections.setSliceLabel("Stdev", n);
            imsProjections.setProcessor(fp, n++);
        }

    }

    // helpers function for radiality calculation

    public void runRadiality(ImageStack ims, ImageStack imsRadiality){

        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        RadialityCL rCL = new RadialityCL(w, h, magnification, ringRadius, sensitivity);

        for (int n = 1; n <= nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fps = rCL.calculateRadiality(ims.getProcessor(n), 0, 0);
            FloatProcessor fpRad = fps[0];

            imsRadiality.addSlice(fpRad);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }
    }

    public void runRadiality(ImageStack ims, ImageStack imsRadiality, ImageStack imsInterpolated){

        int w = ims.getWidth();
        int h = ims.getHeight();
        int nSlices = ims.getSize();

        RadialityCL rCL = new RadialityCL(w, h, magnification, ringRadius, sensitivity);

        for (int n = 1; n <= nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fps = rCL.calculateRadiality(ims.getProcessor(n), 0, 0);
            FloatProcessor fpRad = fps[0];
            FloatProcessor fpInterpolated = fps[1];

            imsRadiality.addSlice(fpRad);
            imsInterpolated.addSlice(fpInterpolated);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return;
            }
        }
    }

}
