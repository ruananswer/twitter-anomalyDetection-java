package com.github.ruananswer.anomaly;

import com.github.ruananswer.testUtility.TestCommon;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Created by ruan on 16-4-7.
 */
public class AnomalyDetectionTest {
    private ArrayList<String> fileStlNames = new ArrayList<String>();
    private ArrayList<String> fileAnomsNames = new ArrayList<String>();
    private static final double EqualPercent = 0.8;

    @BeforeClass
    public void readTestSets() {
        String[] files = TestCommon.readTestSets(" anoms.csv");
        for (String a : files)
            fileAnomsNames.add(a);
        files = TestCommon.readTestSets(" stl.csv");
        for (String a : files)
            fileStlNames.add(a);
    }

    @DataProvider
    public Object[][] DateProviderFromFile() throws IOException, ParseException {
        List<Object[]> data = new ArrayList<Object[]>();
        Collections.sort(fileAnomsNames);
        Collections.sort(fileStlNames);
        if (fileAnomsNames == null || fileStlNames == null) return null;
        Assert.assertEquals(fileAnomsNames.size(), fileStlNames.size());
        String fileStl = "", fileAnoms = "";
        for (int fileIndex = 0; fileIndex < fileAnomsNames.size(); ++fileIndex) {
            fileStl = fileStlNames.get(fileIndex);
            fileAnoms = fileAnomsNames.get(fileIndex);
            String[] lines = TestCommon.getResourceAsString(fileAnoms).split("\n");
            int numAnoms = lines.length;
            int[] anomsIndex = new int[numAnoms - 1];
            double[] anomsScore = new double[numAnoms - 1];
            double maxAnoms = 0.015, alpha = 0.05;
            String[] values = lines[0].split(",");
            if (values.length > 3) {
                maxAnoms = Double.parseDouble(values[3]);
                alpha = Double.parseDouble(values[4]);
            }
            for (int i = 1; i < numAnoms; ++i) {
                values = lines[i].split(",");
                anomsIndex[i - 1] = Integer.parseInt(values[1]);
                anomsScore[i - 1] = Double.parseDouble(values[2]);
            }

            lines = TestCommon.getResourceAsString(fileStl).split("\n");
            int numData = lines.length;
            int seasonality = Integer.parseInt(lines[0].split(",")[0]);
            long[] timestamps = new long[numData - 1];
            double[] series = new double[numData - 1];
            for (int i = 1; i < numData; ++i) {
                values = lines[i].split(",");
                if (values[0].contains("-"))
                    timestamps[i - 1] = (int) new java.text.SimpleDateFormat("yyyy-MM-dd HH:mm:ss").parse(values[0]).getTime();
                else
                    timestamps[i - 1] = Integer.parseInt(values[0]);
                String value = values[1];
                if (value.equals("NaN")) {
                    series[i - 1] = Double.NaN;
                } else {
                    series[i - 1] = Double.valueOf(value);
                }
            }
            data.add(new Object[]{fileAnoms, seasonality, maxAnoms, alpha, timestamps, series, anomsIndex, anomsScore});
        }
        return data.toArray(new Object[][]{});
    }

    @Test(dataProvider = "DateProviderFromFile")
    public void testFunctionAnomalyDetection(String filename, int seasonality, double maxAnoms, double alpha, long[] timestamps,
                                             double[] series, int[] anomsIndex, double[] anomsScore) throws IOException {
        System.out.print(filename);

        DetectAnoms.Config config = new DetectAnoms.Config();
        config.setMaxAnoms(maxAnoms);
        config.setNumObsPerPeriod(seasonality);
        config.setAnomsThreshold(1.15);
        config.setAlpha(alpha);
        DetectAnoms detectAnoms = new DetectAnoms(config);

        long lStartTime = System.nanoTime();
        DetectAnoms.ANOMSResult result = detectAnoms.anomalyDetection(timestamps, series);
        long lEndTime = System.nanoTime();
        long difference = lEndTime - lStartTime; // check different
        System.out.println(" time: " + difference);

        testAnomsSimularity(result, anomsIndex, anomsScore);
    }

    private void testAnomsSimularity(DetectAnoms.ANOMSResult res, int[] anomsIndex_r, double[] score_r) {
        boolean res_test = (res == null && anomsIndex_r.length != 0);
        Assert.assertEquals(res_test, false);
        if (res == null && anomsIndex_r.length == 0)
            return;
        long[] resIndex = res.getAnomsIndex();
        int equalCount = 0;
        for (int i = 0, j = 0; i < resIndex.length && j < anomsIndex_r.length;) {
            if (resIndex[i] + 1 == anomsIndex_r[j]) {
                ++equalCount;
                ++i;
                ++j;
            }
            else if (resIndex[i] + 1 > anomsIndex_r[j])
                ++j;
            else
                ++i;
        }
        if (resIndex.length > 0) {
            boolean equal_test = (double) equalCount / (double) Math.min(anomsIndex_r.length, resIndex.length) > this.EqualPercent;
            Assert.assertEquals(equal_test, true);
        }
    }
}
