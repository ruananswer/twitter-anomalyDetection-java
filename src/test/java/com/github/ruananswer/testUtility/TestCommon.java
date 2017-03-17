package com.github.ruananswer.testUtility;

import org.apache.commons.io.IOUtils;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.StringWriter;

/**
 * Created by ruan on 16-4-18.
 */
public class TestCommon {
   public static String getResourceAsString(String resource) throws IOException {
        InputStream is = ClassLoader.getSystemResourceAsStream(resource);
        StringWriter writer = new StringWriter();
        IOUtils.copy(is, writer);
        return writer.toString();
    }

    public static String[] readTestSets(String regex) {
        File f = new File("./src/test/resources/");
        StlFilenameFilter filter = (new StlFilenameFilter(regex));
        String[] files = f.list(filter);
        return files;
    }

    public static void removeMissingValuesByAveragingNeighbors(double[] arr) {
        for (int i = 0; i < arr.length; i++) {
            if (Double.isNaN(arr[i])) {
                double sum = 0.0;
                int count = 0;
                if (i - 1 >= 0 && !Double.isNaN(arr[i - 1])) {
                    sum += arr[i - 1];
                    count++;
                }
                if (i + 1 < arr.length && !Double.isNaN(arr[i + 1])) {
                    sum += arr[i + 1];
                    count++;
                }
                arr[i] = sum / count;
            }
        }
    }
}
