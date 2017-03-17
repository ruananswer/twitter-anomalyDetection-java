package com.github.ruananswer.stl;

import com.github.ruananswer.testUtility.TestCommon;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.stat.inference.TTest;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by ruan on 16-4-6.
 */

public class STLTest {
    private ArrayList<String> filenames = new ArrayList<String>();
    private static final double SimularityPearsonThreshold = 0.8;

    @BeforeClass
    public void readTestSets() {
        String[] files = TestCommon.readTestSets(" stl.csv");
        for (String a : files)
            filenames.add(a);
    }

    @DataProvider
    public Object[][] DateProviderFromFile() throws IOException {
        List<Object[]> data = new ArrayList<Object[]>();
        for (String filename : filenames) {
            String[] lines = TestCommon.getResourceAsString(filename).split("\n");
            int numData = lines.length;
            int seasonality = Integer.parseInt(lines[0].split(",")[0]);
            long[] timestamps = new long[numData - 1];
            double[] series = new double[numData - 1];
            double[] trend_r = new double[numData - 1];
            double[] residual_r = new double[numData - 1];
            double[] seasonal_r = new double[numData - 1];
            for (int i = 1; i < numData; ++i) {
                String[] values = lines[i].split(",");
                timestamps[i - 1] = Integer.parseInt(values[0]);
                String value = values[1];
                if (value.equals("NA")) {
                    series[i - 1] = Double.NaN;
                } else {
                    series[i - 1] = Double.valueOf(value);
                }
                trend_r[i - 1] = Double.valueOf(values[2]);
                seasonal_r[i - 1] = Double.valueOf(values[3]);
                residual_r[i - 1] = Double.valueOf(values[4]);
            }
            TestCommon.removeMissingValuesByAveragingNeighbors(series);
            numData -= 1;
            data.add(new Object[] {filename, numData, seasonality, timestamps, series, trend_r, residual_r, seasonal_r});
        }
        return data.toArray(new Object[][]{});
    }

    @Test
    public void stlTestNumber() {
        Assert.assertEquals(filenames.size(), 22);
    }

    @Test(dataProvider = "DateProviderFromFile")
    public void testFunctionSTL(String filename, int numData, int seasonality, long[] timestamps, double[] series, double[] trend_r, double[] residual_r, double[] seasonal_r)
    throws IOException {
        System.out.print(filename);

        STLDecomposition.Config config = new STLDecomposition.Config();
        config.setNumberOfDataPoints(numData);
        config.setNumObsPerPeriod(seasonality);
        // if robust
        config.setNumberOfInnerLoopPasses(1);
        config.setNumberOfRobustnessIterations(15);

        STLDecomposition stl = new STLDecomposition(config);

        long lStartTime = System.nanoTime();
        STLResult res = stl.decompose(timestamps, series);
        long lEndTime = System.nanoTime();
        long difference = lEndTime - lStartTime; // check different
        System.out.println(" time: " + difference);

        calSimularity(res, seasonality, trend_r, residual_r, seasonal_r);
    }

    private void calSimularity(STLResult res, int seasonality, double[] trend_r, double[] residual_r, double[] seasonal_r) throws IOException {
        double[] data_trend = res.getTrend();
        double[] data_remainder = res.getRemainder();
        double[] data_seasonal = res.getSeasonal();
        double[] data_seasonal_check = new double[seasonality];
        double[] seasonal_r_check = new double[seasonality];

        double trend_pearson = CheckTimeSeriesSimularityPearson(trend_r, data_trend);
        double remainder_pearson = CheckTimeSeriesSimularityPearson(residual_r, data_remainder);
        double seasonal_pearson = CheckTimeSeriesSimularityPearson(seasonal_r, data_seasonal);
        boolean trend_pearson_test = trend_pearson > SimularityPearsonThreshold;
        boolean remainder_pearson_test = remainder_pearson > SimularityPearsonThreshold;
        boolean seasonal_pearson_test = seasonal_pearson > SimularityPearsonThreshold;
        System.out.println(trend_pearson + "," + remainder_pearson + "," + seasonal_pearson);
        Assert.assertEquals(trend_pearson_test, true);
        Assert.assertEquals(remainder_pearson_test, true);
        Assert.assertEquals(seasonal_pearson_test, true);
        boolean reject_null = checkTimeSeriesSimularityTtest(seasonal_r_check, data_seasonal_check);
        Assert.assertEquals(reject_null, false);
    }

    private boolean checkTimeSeriesSimularityTtest(double[] r_generated, double[] java_generated) {
        TTest t_test_obj = new TTest();
        // Performs a paired t-test evaluating the null hypothesis that the mean of the paired
        // differences between sample1
        // and sample2 is 0 in favor of the two-sided alternative that the mean paired difference is not
        // equal to 0,
        // with significance level 0.05.
        double p_value = t_test_obj.tTest(r_generated, java_generated);
        boolean reject_null_hyphothesis = (p_value < 0.05);
        return reject_null_hyphothesis;
    }

    private double CheckTimeSeriesSimularityPearson(double[] r_generated, double[] java_generated) {
        PearsonsCorrelation pearson = new PearsonsCorrelation();
        double pearson_correlation = pearson.correlation(r_generated, java_generated);
        return pearson_correlation;

    }

}