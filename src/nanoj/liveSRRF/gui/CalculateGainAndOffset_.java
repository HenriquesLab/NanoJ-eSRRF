package nanoj.liveSRRF.gui;

import ij.IJ;
import ij.ImagePlus;
import ij.gui.NonBlockingGenericDialog;
import ij.plugin.PlugIn;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;

public class CalculateGainAndOffset_  implements PlugIn {

    @Override
    public void run(String s) {

        String dirPath = IJ.getDirectory("Choose directory to open...");
        if (dirPath == null) return;

        // need to add extension question here

        ArrayList<String> filesToProcess = this.getFilesToProcess(dirPath, ".tiff");

        for (String fileName: filesToProcess) {

            ImagePlus imp = IJ.openImage(fileName);
            
        }
    }

    private ArrayList<String> getFilesToProcess(String dirPath, String extension) {

        File directory = new File(dirPath);
        if(!directory.exists()){
            return null;
        }

        ArrayList<String> filesToProcess = new ArrayList<String>();

        File[] listOfFilesInFolder = directory.listFiles();
        Arrays.sort(listOfFilesInFolder);

        for (int i=0;i<listOfFilesInFolder.length;i++){
            IJ.showProgress(i+1, listOfFilesInFolder.length);

            if (IJ.escapePressed()) {
                IJ.resetEscape();
                return null;
            }

            if(listOfFilesInFolder[i].isFile() && listOfFilesInFolder[i].toString().endsWith(extension)) {
                String fPath = listOfFilesInFolder[i].getPath();

                filesToProcess.add(fPath);
            }
        }

        return filesToProcess;
    }
}
