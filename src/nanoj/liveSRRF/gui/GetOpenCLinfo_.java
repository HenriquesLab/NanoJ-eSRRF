package nanoj.liveSRRF.gui;

import com.jogamp.common.JogampRuntimeException;
import com.jogamp.opencl.CLDevice;
import com.jogamp.opencl.CLException;
import com.jogamp.opencl.CLPlatform;
import ij.IJ;
import ij.plugin.PlugIn;

public class GetOpenCLinfo_ implements PlugIn {

    public void run(String arg) {

        IJ.log("\\Clear");  // Clear the log window
        IJ.log("-------------------------------------");
        IJ.log("-------------------------------------");

        try{
            CLPlatform.initialize();
        }catch(JogampRuntimeException ex) {
            IJ.log("Could not load Java OpenCL Binding");
            throw new RuntimeException("could not load Java OpenCL Binding");
        }

        CLPlatform[] allPlatforms;
        try{allPlatforms = CLPlatform.listCLPlatforms();}
        catch(CLException ex) {
            IJ.log("Something went wrong initializing OpenCL.");
            throw new RuntimeException("Something went wrong initializing OpenCL.");
        }

        double nFlops;
        double nMaxFlops = 0;
        CLDevice clDeviceMaxFlop = null;

        IJ.log("Number of platform(s): " + allPlatforms.length);
        for (int p = 0; p < allPlatforms.length; p++) {
            IJ.log("-----------------");
            IJ.log("-----------------");
            IJ.log("Platform #"+(p+1)+" ("+allPlatforms[p].getName()+")");
            CLDevice[] allCLdevice = allPlatforms[p].listCLDevices();

            for (int i = 0; i < allCLdevice.length; i++) {
                IJ.log("--------");
                IJ.log("Device #" + i);
                IJ.log("Device name: " + allCLdevice[i].getName());
                IJ.log("Device type: " + allCLdevice[i].getType());
                IJ.log("Max clock: " + allCLdevice[i].getMaxClockFrequency() + " MHz");
                IJ.log("Number of compute units: " + allCLdevice[i].getMaxComputeUnits());

                nFlops = allCLdevice[i].getMaxComputeUnits()*allCLdevice[i].getMaxClockFrequency();
                if (nFlops > nMaxFlops){
                    nMaxFlops = nFlops;
                    clDeviceMaxFlop = allCLdevice[i];
                }
            }
        }

        IJ.log("-----------------");
        assert clDeviceMaxFlop != null;
        IJ.log("Maximum flops device: " + clDeviceMaxFlop.getName());


    }
}
