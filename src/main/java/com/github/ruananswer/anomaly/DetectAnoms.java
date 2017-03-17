package com.github.ruananswer.anomaly;


import com.github.ruananswer.statistics.OnlineNormalStatistics;
import com.github.ruananswer.statistics.QuickMedians;
import com.github.ruananswer.stl.STLDecomposition;
import com.github.ruananswer.stl.STLResult;
import org.apache.commons.math3.distribution.TDistribution;

import java.util.ArrayList;
import java.util.Comparator;

/**
 * Implementation of the underlying algorithm – referred to as Seasonal Hybrid ESD (S-H-ESD) builds upon the Generalized ESD test for detecting anomalies.
 * Created by on 16-4-6.
 */
public class DetectAnoms {
    private final Config config;

    public DetectAnoms(Config config) {
        this.config = config;
    }

    /**
     * The main parameter of function anomalyDetection
     */
    public static class Config {
        /** Maximum number of anomalies that S-H-ESD will detect as a percentage of the data. */
        private double maxAnoms = 0.49;
        /** Defines the number of observations in a single period, and used during seasonal decomposition. */
        private int numObsPerPeriod = 1440;
        /** noms_threshold  use the threshold to filter the anoms. such as if anoms_threshold = 1.05,
         #' then we will filter the anoms that exceed the exceptional critical value 100%-105% */
        private double anomsThreshold = 1.05;
        /** The level of statistical significance with which to accept or reject anomalies. */
        private double alpha = 0.05;

        public Config() {}

        public double getMaxAnoms() {
            return maxAnoms;
        }

        public int getNumObsPerPeriod() {
            return numObsPerPeriod;
        }

        public double getAnomsThreshold() {
            return anomsThreshold;
        }

        public double getAlpha() {
            return alpha;
        }

        public void setMaxAnoms(double maxAnoms) {
            this.maxAnoms = maxAnoms;
        }

        public void setNumObsPerPeriod(int numObsPerPeriod) {
            this.numObsPerPeriod = numObsPerPeriod;
        }

        public void setAnomsThreshold(double anomsThreshold) {
            this.anomsThreshold = anomsThreshold;
        }

        public void setAlpha(double alpha) {
            this.alpha = alpha;
        }
    }

    /**
     * The detected anomalies in a time series using S-H-ESD
     * A list containing the anomalies (anoms) and decomposition components (stl).
     */
    public class ANOMSResult {
        private final long[] anomsIndex;
        private final double[] anomsScore;
        private final double[] dataDecomp;

        private ANOMSResult(long[] anomsIdx, double[] anomsSc, double[] dataDe) {
            this.anomsIndex = anomsIdx;
            this.anomsScore = anomsSc;
            this.dataDecomp = dataDe;
        }

        public long[] getAnomsIndex() {
            return anomsIndex;
        }

        public double[] getAnomsScore() {
            return anomsScore;
        }

        public double[] getDataDecomp() {
            return dataDecomp;
        }
    }

    /** Args:
    #'	series: Time series to perform anomaly detection on.
    #'	max_Anoms: Maximum number of anomalies that S-H-ESD will detect as a percentage of the data.
    #'	alpha: The level of statistical significance with which to accept or reject anomalies.
    #'	num_obs_per_period: Defines the number of observations in a single period, and used during seasonal decomposition.
    #'  Returns:
    #'  A list containing the anomalies (anoms) and decomposition components (stl).
    */
    private ANOMSResult detectAnoms(long[] timestamps, double[] series) {
        if (series == null || series.length < 1) {
            throw new IllegalArgumentException("must supply period length for time series decomposition");
        }
        int numberOfObservations = series.length;
        /**
         * use StlDec function
         * first interpolation to solve problem data to much
         */
        // -- Step 1: Decompose data. This returns a univarite remainder which will be used for anomaly detection. Optionally, we might NOT decompose.
        STLResult data = removeSeasonality(timestamps, series, config.getNumObsPerPeriod());
        double[] data_trend = data.getTrend();
        double[] data_seasonal = data.getSeasonal();

//        Mean mean = new Mean();
//        Variance variance = new Variance();
//        Median median = new Median();

        // Remove the seasonal component, and the median of the data to create the univariate remainder
        double[] dataForSHESD = new double[numberOfObservations];
        double[] dataDecomp = new double[numberOfObservations];

        QuickMedians quickMedian = new QuickMedians(series);
        double medianOfSeries = quickMedian.getMedian();//median.evaluate(series);

        // if the data of has no seasonality, directed use the raw_data into function S-H-ESD !!!
        for (int i = 0; i < numberOfObservations; ++i) {
            dataForSHESD[i] = series[i] - data_seasonal[i] - medianOfSeries;
            dataDecomp[i] = data_trend[i] + data_seasonal[i];
        }
        // Maximum number of outliers that S-H-ESD can detect (e.g. 49% of data)
        int maxOutliers = (int)Math.round(numberOfObservations * config.getMaxAnoms());
        if (maxOutliers == 0)
            throw new IllegalArgumentException("You have " + numberOfObservations + " observations in a period, which is too few. Set a higher value");

        long[] anomsIdx = new long[maxOutliers];
        double[] anomsSc = new double[maxOutliers];
        int numAnoms = 0;

        OnlineNormalStatistics stat = new OnlineNormalStatistics(dataForSHESD);
        QuickMedians quickMedian1 = new QuickMedians(dataForSHESD);
        double dataMean = stat.getMean();//mean.evaluate(dataForSHESD);
        double dataMedian = quickMedian1.getMedian();//median.evaluate(dataForSHESD);
        // use mad replace the variance
        // double dataStd = Math.sqrt(stat.getPopulationVariance());//Math.sqrt(variance.evaluate(dataForSHESD));
        double[] tempDataForMad = new double[numberOfObservations];
        for (int i = 0; i < numberOfObservations; ++i) {
            tempDataForMad[i] = Math.abs(dataForSHESD[i] - dataMedian);
        }
        QuickMedians quickMedian2 = new QuickMedians(tempDataForMad);
        double dataStd = quickMedian2.getMedian();

        if (Math.abs(dataStd) <= 1e-10) {
            //return null;
            throw new IllegalArgumentException("The variance of the series data is zero");
        }

        double[] ares = new double[numberOfObservations];
        for (int i = 0; i < numberOfObservations; ++i) {
            ares[i] = Math.abs(dataForSHESD[i] - dataMedian);
            ares[i] /= dataStd;
        }

        // here use std for the following iterative calculate datastd
        dataStd = Math.sqrt(stat.getPopulationVariance());

        int[] aresOrder = getOrder(ares);
        int medianIndex = numberOfObservations / 2;
        int left = 0, right = numberOfObservations - 1;
        int currentLen = numberOfObservations, tempMaxIdx = 0;
        double R = 0.0, p = 0.0;
        for (int outlierIdx = 1; outlierIdx <= maxOutliers; ++outlierIdx) {
            p = 1.0 - config.getAlpha() / (2 * (numberOfObservations - outlierIdx + 1));
            TDistribution tDistribution = new TDistribution(numberOfObservations - outlierIdx - 1);
            double t = tDistribution.inverseCumulativeProbability(p);
            double lambdaCritical = t * (numberOfObservations - outlierIdx) / Math.sqrt((numberOfObservations
            - outlierIdx - 1 + t * t) * (numberOfObservations - outlierIdx + 1));
            if (left >= right) break;
            if (currentLen < 1) break;

            // remove the largest
            if (Math.abs(dataForSHESD[aresOrder[left]] - dataMedian) > Math.abs(dataForSHESD[aresOrder[right]] - dataMedian)) {
                tempMaxIdx = aresOrder[left];
                ++left;
                ++medianIndex;
            } else {
                tempMaxIdx = aresOrder[right];
                --right;
                --medianIndex;
            }
            // get the R
            R = Math.abs((dataForSHESD[tempMaxIdx] - dataMedian) / dataStd);
            // recalculate the dataMean and dataStd
            dataStd = Math.sqrt(((currentLen - 1) * (dataStd * dataStd + dataMean * dataMean) - dataForSHESD[tempMaxIdx] * dataForSHESD[tempMaxIdx] -
                    ((currentLen - 1) * dataMean - dataForSHESD[tempMaxIdx]) * ((currentLen - 1) * dataMean - dataForSHESD[tempMaxIdx]) /
                            (currentLen - 2)) / (currentLen - 2));
            dataMean = (dataMean * currentLen - dataForSHESD[tempMaxIdx]) / (currentLen - 1);
            dataMedian = dataForSHESD[aresOrder[medianIndex]];
            --currentLen;

            // record the index
            anomsIdx[outlierIdx - 1] = tempMaxIdx;
            anomsSc[outlierIdx - 1] = R;
            if (R < lambdaCritical * config.getAnomsThreshold() || Double.isNaN(dataStd) || Math.abs(dataStd) <= 1e-10) {
                break;
            }
            numAnoms = outlierIdx;
        }
        if (numAnoms > 0) {
            ArrayList<Pair> map = new ArrayList<Pair>();
            for (int i = 0; i < numAnoms; ++i) {
                map.add(new Pair((int)anomsIdx[i], anomsSc[i]));
            }
            map.sort(new PairKeyComparator());
            long[] idx = new long[numAnoms];
            double[] anoms = new double[numAnoms];
            for (int i = 0; i < numAnoms; ++i) {
                idx[i] = map.get(i).key;
                anoms[i] = map.get(i).value;
            }
            return new ANOMSResult(idx, anoms, dataDecomp);
        }
        else
            return null;
    }

    /**
     #' @name AnomalyDetectionTs
     #' @param timestamps & series Time series as a two column data frame where the first column consists of the
     #' timestamps and the second column consists of the observations.
     #' @param max_anoms Maximum number of anomalies that S-H-ESD will detect as a percentage of the
     #' data.
     #' @param numObsPerPeriod the numbers point in one period
     #' @param anoms_threshold  use the threshold to filter the anoms. such as if anoms_threshold = 1.05,
     #' then we will filter the anoms that exceed the exceptional critical value 100%-105%
     #' @param alpha The level of statistical significance with which to accept or reject anomalies.
     */
    public ANOMSResult anomalyDetection(long[] timestamps, double[] series) {
        // Sanity check all input parameters
        if (timestamps == null || timestamps.length < 1 || series == null || series.length < 1 || timestamps.length != series.length)
            throw new IllegalArgumentException("The data is empty or has no equal length.");
        if (config.getMaxAnoms() > 0.49) {
            throw new IllegalArgumentException("max_anoms must be less than 50% of the data points.");
        } else if (config.getMaxAnoms() <= 0) {
            throw new IllegalArgumentException("max_anoms must be positive.");
        }
        /**
         * Main analysis: perform S-H-ESD
         */
        int numberOfObservations = series.length;
        if (config.getMaxAnoms() < 1.0 / numberOfObservations)
            config.setMaxAnoms(1.0 / numberOfObservations);
        removeMissingValuesByAveragingNeighbors(series);

        return detectAnoms(timestamps, series);
    }

    private STLResult removeSeasonality(long[] timestamps, double[] series, int seasonality) {
        STLDecomposition.Config config = new STLDecomposition.Config();
        config.setNumObsPerPeriod(seasonality);
        config.setNumberOfDataPoints(timestamps.length);
        // if robust
        config.setNumberOfInnerLoopPasses(1);
        config.setNumberOfRobustnessIterations(15);

        STLDecomposition stl = new STLDecomposition(config);
        STLResult res = stl.decompose(timestamps, series);

        return res;
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
                if (count != 0)
                    arr[i] = sum / count;
                else
                    arr[i] = 0.0;
            }
        }
    }

    class Pair {
        int key;
        double value;
        public Pair(int k, double v) {
            key = k;
            value = v;
        }
    }

    class PairKeyComparator implements Comparator<Pair> {
        public int compare(Pair a, Pair b) {
            if (a.key != b.key)
                return a.key - b.key;
            return (a.value - b.value) > 0.0 ? 1 : -1;
        }
    }

    class PairValueComparator implements Comparator<Pair> {
        public int compare(Pair a, Pair b) {
            if (a.value != b.value) {
                return (a.value - b.value) > 0 ? -1 : 1;
            }
            return b.key - a.key;
        }
    }

    private int[] getOrder(double[] data) {
        if (data == null || data.length < 1)
            return null;
        int len = data.length;
        ArrayList<Pair> map = new ArrayList<Pair>();
        for (int i = 0; i < len; ++i) {
            map.add(new Pair(i, data[i]));
        }
        map.sort(new PairValueComparator());
        int[] returnOrder = new int[len];
        for (int i = 0; i < len; ++i) {
            returnOrder[i] = map.get(i).key;
        }
        return returnOrder;
    }
}
