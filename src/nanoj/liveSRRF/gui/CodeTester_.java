package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;
import scala.Int;


public class CodeTester_ implements PlugIn {



    public void run(String arg) {

        IJ.log("Float.SIZE: "+(Float.SIZE));
        IJ.log("Integer.SIZE: "+(Integer.SIZE));
        IJ.log("Double.SIZE: "+(Double.SIZE));
        IJ.log("Byte.SIZE: "+(Byte.SIZE));
        IJ.log("Short.SIZE: "+(Short.SIZE));


    }
}
