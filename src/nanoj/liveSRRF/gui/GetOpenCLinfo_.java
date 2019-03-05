package nanoj.liveSRRF.gui;

import com.jogamp.opencl.CLContext;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLPlatform;
import ij.IJ;
import ij.plugin.PlugIn;

public class GetOpenCLinfo_ implements PlugIn {

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        CLPlatform[] allPlatforms = CLPlatform.listCLPlatforms();
        IJ.log("Number of platform(s): " + allPlatforms.length);
        for (int p = 0; p < allPlatforms.length; p++) {
            IJ.log("-----------------");
            IJ.log("-----------------");
            IJ.log("Platform #"+(p+1)+" ("+allPlatforms[p].getProfile()+")");
            //CLContext context = CLContext.create(allPlatforms[p]);

            //CLDevice[] allCLdevice = context.getDevices();
            CLDevice[] allCLdevice = allPlatforms[p].listCLDevices();


            for (int i = 0; i < allCLdevice.length; i++) {
                IJ.log("--------");
                IJ.log("Device #" + i);
                IJ.log("Device name: " + allCLdevice[i].getName());
                IJ.log("Device type: " + allCLdevice[i].getType());
                IJ.log("Max clock: " + allCLdevice[i].getMaxClockFrequency() + " MHz");
                IJ.log("Number of compute units: " + allCLdevice[i].getMaxComputeUnits());
            }

            //context.release();
        }
        IJ.log("-----------------");

    }
}
