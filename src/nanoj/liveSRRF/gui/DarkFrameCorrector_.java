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

        // Sort out dark frames
        if(testMode==choices[0]){
            imsDarkProjections = new ImageStack(wDark, hDark, nDarkProjections);
            FloatProcessor fp;
            int n=1;

            if(averageDark){
                fp = projector.do2DProjection(imsDark, Projections2D.AVERAGE);
                imsDarkProjections.setProcessor(fp, n++);
            }
            if(maximumDark){
                fp = projector.do2DProjection(imsDark, Projections2D.MAX);
                imsDarkProjections.setProcessor(fp, n++);
            }
            if(minimumDark){
                fp = projector.do2DProjection(imsDark, Projections2D.MIN);
                imsDarkProjections.setProcessor(fp, n++);
            }
            if(stdDark){
                fp = projector.do2DProjection(imsDark, Projections2D.STDDEV);
                imsDarkProjections.setProcessor(fp, n++);
            }
        }
        else{

            ImageStack imsDarkRadiality = new ImageStack(wDark*magnification, hDark*magnification);
            ImageStack imsDarkInterpolated = new ImageStack(wDark*magnification, hDark*magnification);

            RadialityCL rCL = new RadialityCL(wDark, hDark, magnification, ringRadius, sensitivity);

            IJ.showStatus("Peforming dark frame radiality calculations...");

            for (int n=1; n<=nSlicesDark; n++) {
                IJ.showProgress(n, nSlicesDark);
                FloatProcessor[] fps = rCL.calculateRadiality(imsDark.getProcessor(n), 0, 0);
                FloatProcessor fpRad = fps[0];
                FloatProcessor fpInterpolated = fps[1];

                imsDarkRadiality.addSlice(fpRad);
                imsDarkInterpolated.addSlice(fpInterpolated);

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return;
                }
             }

            imsDarkProjections = new ImageStack(wDark*magnification, hDark*magnification, nDarkProjections);
            FloatProcessor fp;
            int n=1;

            if(doIW_dark){
                ImageStack imsDarkRadiality_IW = calculateIntensityWeighting(imsDarkRadiality, imsDarkInterpolated);
                if(averageDark){
                    fp = projector.do2DProjection(imsDarkRadiality_IW, Projections2D.AVERAGE);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(maximumDark){
                    fp = projector.do2DProjection(imsDarkRadiality_IW, Projections2D.MAX);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(minimumDark){
                    fp = projector.do2DProjection(imsDarkRadiality_IW, Projections2D.MIN);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(stdDark){
                    fp = projector.do2DProjection(imsDarkRadiality_IW, Projections2D.STDDEV);
                    imsDarkProjections.setProcessor(fp, n++);
                }
            }
            else{
                if(averageDark){
                    fp = projector.do2DProjection(imsDarkRadiality, Projections2D.AVERAGE);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(maximumDark){
                    fp = projector.do2DProjection(imsDarkRadiality, Projections2D.MAX);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(minimumDark){
                    fp = projector.do2DProjection(imsDarkRadiality, Projections2D.MIN);
                    imsDarkProjections.setProcessor(fp, n++);
                }
                if(stdDark){
                    fp = projector.do2DProjection(imsDarkRadiality, Projections2D.STDDEV);
                    imsDarkProjections.setProcessor(fp, n++);
                }
            }

        }

        // end dark frame prep

        ImageStack imsRadialityProjections = new ImageStack(wData*magnification, hData*magnification, nDarkProjections * nRadialityProjections);

        if(testMode==choices[0]) {

            for (int i = 0; i < nDarkProjections; i++) {
                ImageStack imsSubtracted = imsData.duplicate();
                subtract(imsSubtracted, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));

                RadialityCL rCL = new RadialityCL(wData, hData, magnification, ringRadius, sensitivity);

                IJ.showStatus("Peforming data frame radiality calculations for dark projection " + (i + 1) + "/" + nDarkProjections + "...");

                ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
                ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

                for (int n = 1; n <= nSlicesData; n++) {
                    IJ.showProgress(n, nSlicesData);
                    FloatProcessor[] fps = rCL.calculateRadiality(imsSubtracted.getProcessor(n), 0, 0);
                    FloatProcessor fpRad = fps[0];
                    FloatProcessor fpInterpolated = fps[1];

                    imsDataRadiality.addSlice(fpRad);
                    imsDataInterpolated.addSlice(fpInterpolated);

                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        return;
                    }
                }

                FloatProcessor fp;
                int n = i * nRadialityProjections + 1;

                if (doIW_data) {
                    ImageStack imsDataRadiality_IW = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                    if (averageRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality_IW, Projections2D.AVERAGE);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (maximumRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality_IW, Projections2D.MAX);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (minimumRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality_IW, Projections2D.MIN);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (stdRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality_IW, Projections2D.STDDEV);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                } else {
                    if (averageRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality, Projections2D.AVERAGE);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (maximumRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality, Projections2D.MAX);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (minimumRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality, Projections2D.MIN);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                    if (stdRadiality) {
                        fp = projector.do2DProjection(imsDataRadiality, Projections2D.STDDEV);
                        imsRadialityProjections.setProcessor(fp, n++);
                    }
                }

            }
        }

        else if(testMode==choices[1]){

            RadialityCL rCL = new RadialityCL(wData, hData, magnification, ringRadius, sensitivity);

            ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
            ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

            for (int n = 1; n <= nSlicesData; n++) {
                IJ.showProgress(n, nSlicesData);
                FloatProcessor[] fps = rCL.calculateRadiality(imsData.getProcessor(n), 0, 0);
                FloatProcessor fpRad = fps[0];
                FloatProcessor fpInterpolated = fps[1];

                imsDataRadiality.addSlice(fpRad);
                imsDataInterpolated.addSlice(fpInterpolated);

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return;
                }
            }

            ImageStack imsSubtracted = null;

            for (int i = 0; i < nDarkProjections; i++) {

                if(doIW_data) imsSubtracted = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                else imsSubtracted = imsData.duplicate();

                subtract(imsSubtracted, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));

                imsRadialityProjections = new ImageStack(wData * magnification, hData * magnification, nRadialityProjections);
                FloatProcessor fp;
                int n = i * nRadialityProjections + 1;

                if (averageRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.AVERAGE);
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (maximumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MAX);
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (minimumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MIN);
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (stdRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.STDDEV);
                    imsRadialityProjections.setProcessor(fp, n++);
                }

            }

        }

        else if(testMode==choices[2]){

            RadialityCL rCL = new RadialityCL(wData, hData, magnification, ringRadius, sensitivity);

            ImageStack imsDataRadiality = new ImageStack(wDark * magnification, hDark * magnification);
            ImageStack imsDataInterpolated = new ImageStack(wDark * magnification, hDark * magnification);

            for (int n = 1; n <= nSlicesData; n++) {
                IJ.showProgress(n, nSlicesData);
                FloatProcessor[] fps = rCL.calculateRadiality(imsData.getProcessor(n), 0, 0);
                FloatProcessor fpRad = fps[0];
                FloatProcessor fpInterpolated = fps[1];

                imsDataRadiality.addSlice(fpRad);
                imsDataInterpolated.addSlice(fpInterpolated);

                if (IJ.escapePressed()) {
                    IJ.resetEscape();
                    return;
                }
            }

            ImageStack imsSubtracted = null;

            for(int i=0; i<nDarkProjections; i++){
                if(doIW_data) imsSubtracted = calculateIntensityWeighting(imsDataRadiality, imsDataInterpolated);
                else imsSubtracted = imsData.duplicate();

                FloatProcessor fp;
                int n = i * nRadialityProjections + 1;

                if (averageRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.AVERAGE);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (maximumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MAX);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (minimumRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.MIN);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
                if (stdRadiality) {
                    fp = projector.do2DProjection(imsSubtracted, Projections2D.STDDEV);
                    subtract(fp, (FloatProcessor) imsDarkProjections.getProcessor(i + 1));
                    imsRadialityProjections.setProcessor(fp, n++);
                }
            }

        }

        new ImagePlus("Subtracted stacks - "+testMode, imsRadialityProjections).show();


    }


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

        ImageStack imsSubtracted = new ImageStack(w, h, nSlices);

        for(int n=1; n<=nSlices; n++){
            float[] imsPixels = (float[]) ims.getProcessor(n).convertToFloatProcessor().getPixels();
            float[] fpPixels = (float[]) fp.getPixels();
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

}
