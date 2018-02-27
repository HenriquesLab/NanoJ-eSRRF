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

    String testMode, titleDark, titleData;
    boolean doIW_dark, doIW_data;
    boolean averageDark, maximumDark, minimumDark, stdDark;
    int magnification;
    double ringRadius, sensitivity;
    boolean averageRadiality, maximumRadiality, minimumRadiality, stdRadiality;

    ImagePlus impDark, impData;
    Projections2D projector = new Projections2D();

    @Override
    public void run(String arg) {

        String[] choices = new String[]{"Case 1: Raw data subtraction", "Case 2: Pre-projection subtraction",
                                        "Case 3: Post-projection subtraction"};
        String[] imageTitles = WindowManager.getImageTitles();

        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Dark frame correction");

        gd.addChoice("Test mode", choices, choices[0]);

        gd.addMessage("Dark frames settings");
        gd.addChoice("Dark frames image", imageTitles, imageTitles[0]);
        gd.addCheckbox("Do intensity weighting", getPrefs("doIW_dark", true));

        gd.addMessage("Dark frames projection options");
        gd.addCheckbox("Average", getPrefs("averageDark", true));
        gd.addCheckbox("Maximum Intensity", getPrefs("maximumDark", true));
        gd.addCheckbox("Minimum Intensity", getPrefs("minimumDark", true));
        gd.addCheckbox("Standard Deviation", getPrefs("stdDark", true));

        gd.addMessage("Data frames settings");
        gd.addChoice("Data frames image", imageTitles, imageTitles[0]);
        gd.addCheckbox("Do intensity weighting", getPrefs("doIW_data", true));

        gd.addMessage("Radiality settings");
        gd.addNumericField("Magnification", getPrefs("magnification", 5),0);
        gd.addNumericField("Ring radius", getPrefs("ringRadius", 0.5f), 1);
        gd.addNumericField("Radiality sensitivity", getPrefs("sensitivity", 2), 0);

        gd.addMessage("Radiality projection options");
        gd.addCheckbox("Average", getPrefs("averageRadiality", true));
        gd.addCheckbox("Maximum Intensity", getPrefs("maximumRadiality", true));
        gd.addCheckbox("Minimum Intensity", getPrefs("minimumRadiality", true));
        gd.addCheckbox("Standard Deviation", getPrefs("stdRadiality", true));

        gd.showDialog();
        if (gd.wasCanceled()) return;

        testMode = gd.getNextChoice();

        titleDark = gd.getNextChoice();

        doIW_dark = gd.getNextBoolean();
        setPrefs("doIW_dark", doIW_dark);

        averageDark = gd.getNextBoolean();
        setPrefs("averageDark", averageDark);
        maximumDark = gd.getNextBoolean();
        setPrefs("maxmimumDark", maximumDark);
        minimumDark = gd.getNextBoolean();
        setPrefs("minimumDark", minimumDark);
        stdDark = gd.getNextBoolean();
        setPrefs("stdDark", stdDark);

        titleData = gd.getNextChoice();

        doIW_data = gd.getNextBoolean();
        setPrefs("doIW_data", doIW_data);

        magnification = (int) gd.getNextNumber();
        setPrefs("magnification", (float) magnification);
        ringRadius = gd.getNextNumber();
        setPrefs("ringRadius", (float) ringRadius);
        sensitivity = gd.getNextNumber();
        setPrefs("sensitivity", (float) sensitivity);

        averageRadiality = gd.getNextBoolean();
        setPrefs("averageRadiality", averageRadiality);
        maximumRadiality = gd.getNextBoolean();
        setPrefs("maxmimumRadiality", maximumRadiality);
        minimumRadiality = gd.getNextBoolean();
        setPrefs("minimumRadiality", minimumRadiality);
        stdRadiality = gd.getNextBoolean();
        setPrefs("stdRadiality", stdRadiality);

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

        assert (wDark==wData && hDark==hData);

        boolean[] darkProjections = new boolean[]{averageDark, maximumDark, minimumDark, stdDark};
        boolean[] radialityProjections = new boolean[]{averageRadiality, maximumRadiality, minimumRadiality,
                                                        stdRadiality};

        int nDarkProjections = 0;
        for(int i=0; i<darkProjections.length; i++){
            if(darkProjections[i]) nDarkProjections++;
        }

        int nRadialityProjections = 0;
        for(int i=0; i<radialityProjections.length; i++){
            if(radialityProjections[i]) nRadialityProjections++;
        }

        ImageStack imsDarkProjections;
        String[] darkLabels = new String[nDarkProjections];

        // Sort out dark frames
        if(testMode==choices[0]){
            imsDarkProjections = new ImageStack(wDark, hDark, nDarkProjections);
            doZProjections(imsDarkProjections, imsDark, 1, false, new String[0]);
            for(int n=0; n<nDarkProjections; n++){
                darkLabels[n] = imsDarkProjections.getSliceLabel(n+1);
            }
        }
        else{

            ImageStack imsDarkRadiality = new ImageStack(wDark*magnification, hDark*magnification);
            ImageStack imsDarkInterpolated = new ImageStack(wDark*magnification, hDark*magnification);

            IJ.showStatus("Performing dark frame radiality calculations...");
            runRadiality(imsDark, imsDarkRadiality, imsDarkInterpolated);

            imsDarkProjections = new ImageStack(wDark*magnification, hDark*magnification, nDarkProjections);

            if(doIW_dark){
                ImageStack imsDarkRadiality_IW = calculateIntensityWeighting(imsDarkRadiality, imsDarkInterpolated);
                doZProjections(imsDarkProjections, imsDarkRadiality_IW, 1, false, new String[0]);
            }
            else{
                doZProjections(imsDarkProjections, imsDarkRadiality, 1, false, new String[0]);
            }

            for(int n=0; n<nDarkProjections; n++){
                darkLabels[n] = imsDarkProjections.getSliceLabel(n+1);
            }

        }

        // end dark frame prep

        ImageStack imsRadialityProjections = new ImageStack(wData*magnification, hData*magnification, nDarkProjections * nRadialityProjections);

        if(testMode==choices[0]) {

            for (int i = 0; i < nDarkProjections; i++) {

                ImageStack imsSubtracted = imsData.duplicate();
                subtract(imsSubtracted, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));

                ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
                ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

                runRadiality(imsSubtracted, imsDataRadiality, imsDataInterpolated);

                IJ.showStatus("Performing data frame radiality calculations for dark projection " + (i + 1) + "/" + nDarkProjections + "...");

                int n = i * nRadialityProjections + 1;

                if (doIW_data) {
                    ImageStack imsDataRadiality_IW = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                    doZProjections(imsRadialityProjections, imsDataRadiality_IW, n, true, darkLabels);
                }
                else {
                    doZProjections(imsRadialityProjections, imsDataRadiality, n, true, darkLabels);
                }

            }
        }

        else if(testMode==choices[1]){

            ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
            ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

            runRadiality(imsData, imsDataRadiality, imsDataInterpolated);

            ImageStack imsSubtracted;

            for (int i = 0; i < nDarkProjections; i++) {

                if(doIW_data){
                    imsSubtracted = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                }
                else imsSubtracted = imsData.duplicate();

                subtract(imsSubtracted, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));

                int n = i * nRadialityProjections + 1;
                doZProjections(imsRadialityProjections, imsSubtracted, n, true, darkLabels);

            }

        }

        else if(testMode==choices[2]){

            ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
            ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

            runRadiality(imsData, imsDataRadiality, imsDataInterpolated);

            ImageStack imsSubtracted;

            for(int i=0; i<nDarkProjections; i++){

                String darkSliceLabel = imsDarkProjections.getSliceLabel(i+1);
                String IWLabel = "";

                if(doIW_data){
                    imsSubtracted = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                    IWLabel = "IW-";
                }
                else imsSubtracted = imsData.duplicate();

                FloatProcessor fp;
                int n = i * nRadialityProjections + 1;

                if (averageRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.AVERAGE);
                    imsRadialityProjections.setSliceLabel(darkSliceLabel+"_RADIALITY-"+IWLabel+"AVERAGE", n);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (maximumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MAX);
                    imsRadialityProjections.setSliceLabel(darkSliceLabel+"_RADIALITY-"+IWLabel+"MAX", n);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (minimumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MIN);
                    imsRadialityProjections.setSliceLabel(darkSliceLabel+"_RADIALITY-"+IWLabel+"MIN", n);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (stdRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.STDDEV);
                    imsRadialityProjections.setSliceLabel(darkSliceLabel+"_RADIALITY-"+IWLabel+"STDEV", n);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
            }

        }

        new ImagePlus("Subtracted stacks - "+testMode, imsRadialityProjections).show();


    }


    //helper functions for prefs handling

    public float getPrefs(String key, float defaultValue) {
        return (float) prefs.get(prefsHeader+"."+key, defaultValue);
    }

    public boolean getPrefs(String key, boolean defaultValue) {
        return prefs.get(prefsHeader+"."+key, defaultValue);
    }

    public void setPrefs(String key, float defaultValue) {
        prefs.set(prefsHeader+"."+key, defaultValue);
    }

    public void setPrefs(String key, boolean defaultValue) {
        prefs.set(prefsHeader+"."+key, defaultValue);
    }

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

    public void subtract(FloatProcessor fp1, FloatProcessor fp2){
        int w = fp1.getWidth();
        int h = fp1.getHeight();

        float[] fp1Pixels = (float[]) fp1.getPixels();
        float[] fp2Pixels = (float[]) fp2.getPixels();
        for(int i=0; i<w*h; i++){
            fp1Pixels[i] -= fp2Pixels[i];
        }
    }

    // helper functions for z-projecting

    public void doZProjections(ImageStack imsProjections, ImageStack imsRaw, int startIndex, boolean doingRadiality,
                               String[] darkLabels){

        int n=startIndex;
        FloatProcessor fp;
        boolean doAverage, doMaximum, doMinimum, doStd;
        String label;

        if(doingRadiality){
            doAverage = averageRadiality;
            doMaximum = maximumRadiality;
            doMinimum = minimumRadiality;
            doStd = stdRadiality;
        }
        else{
            doAverage = averageDark;
            doMaximum = maximumDark;
            doMinimum = minimumDark;
            doStd = stdDark;
        }

        int labelIndex = 0;

        if (doAverage) {
            fp = projector.do2DProjection(imsRaw, Projections2D.AVERAGE);
            if(doingRadiality){
                label = makeRadialitySliceLabel(Projections2D.AVERAGE, darkLabels[labelIndex++]);
            }
            else{
                label = makeDarkSliceLabel(Projections2D.AVERAGE);
            }
            imsProjections.setSliceLabel(label, n);
            imsProjections.setProcessor(fp, n++);
        }
        if (doMaximum) {
            fp = projector.do2DProjection(imsRaw, Projections2D.MAX);
            if(doingRadiality){
                label = makeRadialitySliceLabel(Projections2D.MAX, darkLabels[labelIndex++]);
            }
            else{
                label = makeDarkSliceLabel(Projections2D.MAX);
            }
            imsProjections.setSliceLabel(label, n);
            imsProjections.setProcessor(fp, n++);
        }
        if (doMinimum) {
            fp = projector.do2DProjection(imsRaw, Projections2D.MIN);
            if(doingRadiality){
                label = makeRadialitySliceLabel(Projections2D.MIN, darkLabels[labelIndex++]);
            }
            else{
                label = makeDarkSliceLabel(Projections2D.MIN);
            }
            imsProjections.setSliceLabel(label, n);
            imsProjections.setProcessor(fp, n++);
        }
        if (doStd) {
            fp = projector.do2DProjection(imsRaw, Projections2D.STDDEV);
            if(doingRadiality){
                label = makeRadialitySliceLabel(Projections2D.STDDEV, darkLabels[labelIndex++]);
            }
            else{
                label = makeDarkSliceLabel(Projections2D.STDDEV);
            }
            imsProjections.setSliceLabel(label, n);
            imsProjections.setProcessor(fp, n++);
        }

    }

    public String makeDarkSliceLabel(int projectionType){
        String projectionString, IWString = "";

        if(projectionType==Projections2D.AVERAGE) projectionString = "AVG";
        else if (projectionType==Projections2D.MAX) projectionString = "MAX";
        else if (projectionType==Projections2D.MIN) projectionString = "MIN";
        else projectionString = "STD";

        if(doIW_dark) IWString = "IW-";

        return "Dark-"+IWString+projectionString;
    }

    public String makeRadialitySliceLabel(int projectionType, String darkLabel){
        String projectionString, IWString = "";

        if(projectionType==Projections2D.AVERAGE) projectionString = "AVG";
        else if (projectionType==Projections2D.MAX) projectionString = "MAX";
        else if (projectionType==Projections2D.MIN) projectionString = "MIN";
        else projectionString = "STD";

        if(doIW_data) IWString = "IW-";

        return darkLabel+"_Radiality-"+IWString+projectionString;
    }

    // helper function for radiality calculation

    public void runRadiality(ImageStack imsRaw, ImageStack imsRadiality, ImageStack imsInterpolated){

        int w = imsRaw.getWidth();
        int h = imsRaw.getHeight();
        int nSlices = imsRaw.getSize();

        RadialityCL rCL = new RadialityCL(w, h, magnification, ringRadius, sensitivity);

        for (int n = 1; n <= nSlices; n++) {
            IJ.showProgress(n, nSlices);
            FloatProcessor[] fps = rCL.calculateRadiality(imsRaw.getProcessor(n), 0, 0);
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
