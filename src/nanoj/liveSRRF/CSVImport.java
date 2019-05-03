package nanoj.liveSRRF;

import ij.IJ;
import org.apache.commons.io.input.CountingInputStream;
import org.scijava.table.DefaultFloatTable;

import java.io.*;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.Locale;
import java.util.Scanner;
import java.io.IOException;

// A large part of this was taken from ThunderSTORM code on Github on 2019-05-01
// https://github.com/zitmen/thunderstorm/blob/master/src/main/java/cz/cuni/lf1/lge/ThunderSTORM/ImportExport/DLMImportExport.java


public class CSVImport {

    public DefaultFloatTable dataTable;
    private String separator = ",";
    private int floatPrecision = 5;

    public CSVImport(String filePath) throws IOException {

        File file = new File(filePath);
        long fileSize = file.length();

        //im using file size and counting bytes read to track progress (without having to know line count in advance)
        Scanner sc = null;
        CountingInputStream countingInputStream = null;


        try {
            countingInputStream = new CountingInputStream(new BufferedInputStream(new FileInputStream(file)));
            sc = new Scanner(countingInputStream, "UTF-8");

            if (!sc.hasNextLine()) return;
            String[] firstLine = splitLine(sc.nextLine());
            if (firstLine.length < 2) return;

            int nCols = firstLine.length;

            dataTable = new DefaultFloatTable(nCols,0);

            DecimalFormat df = new DecimalFormat();
            DecimalFormatSymbols symbols = DecimalFormatSymbols.getInstance(Locale.US);
            symbols.setInfinity("Infinity");
            symbols.setNaN("NaN");
            df.setDecimalFormatSymbols(symbols);
            df.setGroupingUsed(false);
            df.setRoundingMode(RoundingMode.HALF_EVEN);
            df.setMaximumFractionDigits(floatPrecision);

            int nRows = -1;
            float[] values = new float[nCols];
            String[] line;
            while (sc.hasNextLine()) {
                line = splitLine(sc.nextLine());
                if (line.length == 1 && line[0].isEmpty()) continue; // just an empty line...do not report
                if (line.length != nCols) { // but report when there is a corrupted row
                    IJ.log("\\Update:A line has different number of items than the header! Skipping over...");
                    continue;
                }
                nRows++;

                try {
                    for (int c = 0, ci = 0, cm = line.length; c < cm; c++) {
                        values[ci] = df.parse(line[c]).floatValue();
                        ci++;
                    }

                    // insert data in table
                    dataTable.appendRow();
                    for (int i = 0; i<nCols; i++){
                        dataTable.set(i, nRows, values[i]);
                    }

                } catch (ParseException ex) {
                    IJ.log("\\Update:Invalid number format! Skipping over...");
                }

                IJ.showProgress((double) countingInputStream.getByteCount() / (double)fileSize);
            }
        } finally {
            if (countingInputStream != null) countingInputStream.close();
            if (sc != null) sc.close();
        }

    }

    private String[] splitLine(String line) {
        String[] arr = line.split(separator);
        for (int i = 0; i < arr.length; i++) {
            arr[i] = arr[i].trim();
            if (arr[i].startsWith("\"") && arr[i].endsWith("\"") ||
                    arr[i].startsWith("'") && arr[i].endsWith("'") ||
                    arr[i].startsWith("`") && arr[i].endsWith("`")) {
                arr[i] = arr[i].substring(1, arr[i].length() - 1);
            }
        }
        return arr;
    }
}
