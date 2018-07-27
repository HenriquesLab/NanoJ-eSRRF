package nanoj.liveSRRF.gui;

import com.jogamp.opencl.CLDevice;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import nanoj.core2.NanoJPrefs;
import nanoj.liveSRRF.liveSRRF_CL;

import java.awt.*;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class ParametersSweep_ implements PlugIn {

    // Basic formats
    private int magnification,
            nSlices,
            width,
            height,
            blockSize;

    private boolean calculateAVG,
            calculateSTD;

    private float[] fwhmArray;
    private int[] sensitivityArray,
            nframeArray;

    // Image formats
    private ImagePlus imp;

    // Advanced formats
    private NanoJPrefs prefs = new NanoJPrefs(this.getClass().getName());
    liveSRRF_CL liveSRRF;


    public void run(String arg) {


        // Get raw data
        imp = WindowManager.getCurrentImage();
        if (imp == null) imp = IJ.openImage();
        if (imp == null) return;
        imp.show();

        nSlices = imp.getImageStack().getSize();
        width = imp.getImageStack().getWidth();
        height = imp.getImageStack().getHeight();

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");
        IJ.log("liveSRRF - Parameters sweep");
        LocalDateTime now = LocalDateTime.now();
        DateTimeFormatter formatter = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");
        IJ.log(now.format(formatter));
        liveSRRF = new liveSRRF_CL();
        CLDevice[] allDevices = liveSRRF.checkDevices();

        // Initilizaing string for device choice
        String[] deviceNames = new String[allDevices.length + 1];
        deviceNames[0] = "Default device";

        for (int i = 1; i <= allDevices.length; i++) {
            deviceNames[i] = allDevices[i - 1].getName();
        }

        // Build GUI
        Font headerFont = new Font("Arial", Font.BOLD, 16);
        NonBlockingGenericDialog gd = new NonBlockingGenericDialog("liveSRRF - Parameters sweep");
        gd.addMessage("-=-= Fixed SRRF parameters =-=-\n", headerFont);
        gd.addNumericField("Magnification (default: 5)", prefs.get("magnification", 5), 0);

        gd.addMessage("-=-= Sweeping SRRF parameters =-=-\n", headerFont);
        gd.addMessage("FWHM\n", headerFont);
        gd.addNumericField("Start", prefs.get("fwhm0", 2), 2);
        gd.addNumericField("Delta", prefs.get("deltafwhm", 0.5f), 2);
        gd.addNumericField("Number", prefs.get("n_fwhm", 5), 0);

        gd.addMessage("Sensitivity\n", headerFont);
        gd.addNumericField("Start", prefs.get("S0", 1), 0);
        gd.addNumericField("Delta", prefs.get("deltaS", 1), 0);
        gd.addNumericField("Number", prefs.get("n_S", 4), 0);

        gd.addMessage("# frames for SRRF\n", headerFont);
        gd.addNumericField("Start", prefs.get("nf0", 50), 0);
        gd.addNumericField("Delta", prefs.get("deltanf", 25), 0);
        gd.addNumericField("Number", prefs.get("n_nf", 3), 0);

        gd.addMessage("-=-= Reconstructions =-=-\n", headerFont);
        gd.addCheckbox("AVG reconstruction (default: on)", prefs.get("calculateAVG", true));
        gd.addCheckbox("STD reconstruction (default: off)", prefs.get("calculateSTD", false));

        gd.addMessage("-=-= GPU processing =-=-\n", headerFont);
        gd.addNumericField("Analysis block size (default: 20000)", prefs.get("blockSize", 20000), 0);
        gd.addMessage("A large analysis block size will speed up the analysis but will use\n" +
                "more resources and may slow down your computer.");

        gd.showDialog();

        // If the GUI was cancelled
        if (gd.wasCanceled()) {
            liveSRRF.release();
            return;
        }

        grabSettings(gd);

        ImageStack imsBuffer;
        ImageStack imsInt = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFavg = new ImageStack(width * magnification, height * magnification);
        ImageStack imsSRRFstd = new ImageStack(width * magnification, height * magnification);

        ImagePlus impTemp = new ImagePlus();
        impTemp.copyScale(imp); // make sure we copy the pixel sizes correctly across
        Calibration cal = impTemp.getCalibration();
        cal.pixelWidth /= magnification;
        cal.pixelHeight /= magnification;
        cal.setUnit(imp.getCalibration().getUnit());

        ImageStack imsRawData;

        int n_calculation = nframeArray.length * sensitivityArray.length * fwhmArray.length;
        IJ.log("Number of calculations planned: " + n_calculation);
        int r = 1;


        for (int thisnf : nframeArray) {
            for (int thisSensitivity : sensitivityArray) {
                for (float thisfwhm : fwhmArray) {
                    IJ.log("--------");
                    IJ.log("SRRF frame: " + r +"/"+n_calculation);
                    IJ.showProgress(r, n_calculation);

                    // Check if user is cancelling calculation
                    if (IJ.escapePressed()) {
                        IJ.resetEscape();
                        liveSRRF.release();
                        return;
                    }

                    IJ.log("Number of frame for SRRF: "+thisnf);
                    IJ.log("FWHM: "+thisfwhm + " pixels");
                    IJ.log("Sensitivity: "+thisSensitivity);

                    imsRawData = new ImageStack(width, height);
                    for (int f = 1; f <= thisnf; f++) {
                        imp.setSlice(f);
                        imsRawData.addSlice(imp.getProcessor());
                    }

                    liveSRRF.initialise(width, height, magnification, thisfwhm, thisSensitivity, thisnf, thisnf, blockSize,null);
                    liveSRRF.resetFramePosition();
                    liveSRRF.calculateSRRF(imsRawData);
                    imsBuffer = liveSRRF.readSRRFbuffer();

                    if (calculateAVG) imsSRRFavg.addSlice(imsBuffer.getProcessor(1));
                    if (calculateSTD) imsSRRFstd.addSlice(imsBuffer.getProcessor(2));
                    imsInt.addSlice(imsBuffer.getProcessor(3));

                    r++;
                }
            }
        }

        // Release the GPU
        liveSRRF.release();

        //Display results
        if (calculateAVG) {
            ImagePlus impSRRFavg = new ImagePlus(imp.getTitle() + " - liveSRRF (AVG)", imsSRRFavg);
            impSRRFavg.setCalibration(cal);
            IJ.run(impSRRFavg, "Enhance Contrast", "saturated=0.5");
            impSRRFavg.show();
        }

        if (calculateSTD) {
            ImagePlus impSRRFstd = new ImagePlus(imp.getTitle() + " - liveSRRF (STD)", imsSRRFstd);
            impSRRFstd.setCalibration(cal);
            IJ.run(impSRRFstd, "Enhance Contrast", "saturated=0.5");
            impSRRFstd.show();
        }


        ImagePlus impInt = new ImagePlus(imp.getTitle() + " - Interpolated image", imsInt);
        impInt.setCalibration(cal);
        IJ.run(impInt, "Enhance Contrast", "saturated=0.5");
        impInt.show();

        IJ.log("-------------------------------------");
        IJ.log("RAM used: " + IJ.freeMemory());
        IJ.log("Bye-bye !");

    }

    //    -------------------------------------------------------------------------------------
    //    -------------------------- Here lies dragons and functions --------------------------
    //    -------------------------------------------------------------------------------------


    //    --- Grab settings ---
    private void grabSettings(GenericDialog gd) {

        magnification = (int) gd.getNextNumber();
        float fwhm0 = (float) gd.getNextNumber();
        float deltafwhm = (float) gd.getNextNumber();
        int n_fwhm = (int) gd.getNextNumber();

        int S0 = (int) gd.getNextNumber();
        int deltaS = (int) gd.getNextNumber();
        int n_S = (int) gd.getNextNumber();

        int nf0 = (int) gd.getNextNumber();
        int deltanf = (int) gd.getNextNumber();
        int n_nf = (int) gd.getNextNumber();

        calculateAVG = gd.getNextBoolean();
        calculateSTD = gd.getNextBoolean();
        blockSize = (int) gd.getNextNumber();

        fwhmArray = new float[n_fwhm];
        for (int i = 0; i < n_fwhm; i++) {
            fwhmArray[i] = fwhm0 + i * deltafwhm;
        }

        sensitivityArray = new int[n_S];
        for (int i = 0; i < n_S; i++) {
            sensitivityArray[i] = S0 + i * deltaS;
        }

        nframeArray = new int[n_nf];
        for (int i = 0; i < n_nf; i++) {
            nframeArray[i] = nf0 + i * deltanf;
        }


        prefs.set("magnification", (float) magnification);
        prefs.set("fwhm0", fwhm0);
        prefs.set("deltafwhm", deltafwhm);
        prefs.set("n_fwhm", n_fwhm);

        prefs.set("S0", S0);
        prefs.set("deltaS", deltaS);
        prefs.set("n_S", n_S);

        prefs.set("nf0", nf0);
        prefs.set("deltanf", deltanf);
        prefs.set("n_nf", n_nf);

        prefs.set("calculateAVG", calculateAVG);
        prefs.set("calculateSTD", calculateSTD);

        prefs.set("blockSize", blockSize);


        prefs.save();


    }


}
