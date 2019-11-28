package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.plugin.PlugIn;


public class CodeTester_ implements PlugIn {



    public void run(String arg) {

        byte b = 0;
        float f = 0;

        for (int i = 0; i < 20; i++) {
            b += 10;
            f = (float) b;
            IJ.log("b: "+b);
            IJ.log("f: "+f);

        }


    }
}
