package com.github.ruananswer.stl;

import java.util.ArrayList;

/**
 * Implementation of STL: A Seasonal-Trend Decomposition Procedure based on Loess.
 * <p>
 * Robert B. Cleveland et al., "STL: A Seasonal-Trend Decomposition Procedure based on Loess,"
 * in Journal of Official Statistics Vol. 6 No. 1, 1990, pp. 3-73
 * </p>
 * <p>
 * Hafen, R. P. "Local regression models: Advancements, applications, and new methods." (2010).
 * </p>
 * Created by on 16-4-14.
 */
public class STLDecomposition {
    private final Config config;

    public STLDecomposition(Config config) {
        config.check();
        this.config = config;
    }

    public static class Config {
        /** The number of observations in each cycle of the seasonal component, n_p */
        private int numObsPerPeriod = -1;
        /** s.window either the character string \code{"periodic"} or the span (in lags) of the loess window for seasonal extraction,
         * which should be odd.  This has no default. */
        private int sWindow = -1;
        /** s.degree degree of locally-fitted polynomial in seasonal extraction.  Should be 0, 1, or 2. */
        private int sDegree = 1;
        /** t.window the span (in lags) of the loess window for trend extraction, which should be odd.
         *  If \code{NULL}, the default, \code{nextodd(ceiling((1.5*period) / (1-(1.5/s.window))))}, is taken.*/
        private int tWindow = -1;
        /** t.degree degree of locally-fitted polynomial in trend extraction.  Should be 0, 1, or 2. */
        private int tDegree = 1;
        /** l.window the span (in lags) of the loess window of the low-pass filter used for each subseries.
         * Defaults to the smallest odd integer greater than or equal to \code{n.p}
         * which is recommended since it prevents competition between the trend and seasonal components.
         * If not an odd integer its given value is increased to the next odd one.*/
        private int lWindow = -1;
        /** l.degree degree of locally-fitted polynomial for the subseries low-pass filter.  Should be 0, 1, or 2. */
        private int lDegree = 1;
        /** s.jump s.jump,t.jump,l.jump,fc.jump integers at least one to increase speed of the respective smoother.
         * Linear interpolation happens between every \code{*.jump}th value. */
        private int sJump = -1;
        /** t.jump */
        private int tJump = -1;
        /** l.jump */
        private int lJump = -1;
        /** critfreq the critical frequency to use for automatic calculation of smoothing windows for the trend and high-pass filter. */
        private double critFreq = 0.05;
        /** The number of passes through the inner loop, n_i */
        private int numberOfInnerLoopPasses = 2;
        /** The number of robustness iterations of the outer loop, n_o */
        private int numberOfRobustnessIterations = 1;
        /** sub.labels optional vector of length n.p that contains the labels of the subseries in their natural order (such as month name, day of week, etc.),
         * used for strip labels when plotting.  All entries must be unique. */
        private int[] subLabels = null;
        /** the number of series to decompose */
        private int numberOfDataPoints = -1;

        public int getNumObsPerPeriod() {
            return numObsPerPeriod;
        }

        public void setNumObsPerPeriod(int numObsPerPeriod) {
            this.numObsPerPeriod = numObsPerPeriod;
        }

        public int getsWindow() {
            return sWindow;
        }

        public void setsWindow(int sWindow) {
            this.sWindow = sWindow;
        }

        public int getsDegree() {
            return sDegree;
        }

        public void setsDegree(int sDegree) {
            this.sDegree = sDegree;
        }

        public int gettWindow() {
            return tWindow;
        }

        public void settWindow(int tWindow) {
            this.tWindow = tWindow;
        }

        public int gettDegree() {
            return tDegree;
        }

        public void settDegree(int tDegree) {
            this.tDegree = tDegree;
        }

        public int getlWindow() {
            return lWindow;
        }

        public void setlWindow(int lWindow) {
            this.lWindow = lWindow;
        }

        public int getlDegree() {
            return lDegree;
        }

        public void setlDegree(int lDegree) {
            this.lDegree = lDegree;
        }

        public int getsJump() {
            return sJump;
        }

        public void setsJump(int sJump) {
            this.sJump = sJump;
        }

        public int gettJump() {
            return tJump;
        }

        public void settJump(int tJump) {
            this.tJump = tJump;
        }

        public int getlJump() {
            return lJump;
        }

        public void setlJump(int lJump) {
            this.lJump = lJump;
        }

        public double getCritFreq() {
            return critFreq;
        }

        public void setCritFreq(double critFreq) {
            this.critFreq = critFreq;
        }

        public int getNumberOfInnerLoopPasses() {
            return numberOfInnerLoopPasses;
        }

        public void setNumberOfInnerLoopPasses(int numberOfInnerLoopPasses) {
            this.numberOfInnerLoopPasses = numberOfInnerLoopPasses;
        }

        public int getNumberOfRobustnessIterations() {
            return numberOfRobustnessIterations;
        }

        public void setNumberOfRobustnessIterations(int numberOfRobustnessIterations) {
            this.numberOfRobustnessIterations = numberOfRobustnessIterations;
        }

        public int[] getSubLabels() {
            return subLabels;
        }

        public void setSubLabels(int[] subLabels) {
            this.subLabels = subLabels;
        }

        public int getNumberOfDataPoints() {
            return numberOfDataPoints;
        }

        public void setNumberOfDataPoints(int numberOfDataPoints) {
            this.numberOfDataPoints = numberOfDataPoints;
        }

        public Config() {
        }

        public void check() {
            checkPeriodicity(numObsPerPeriod, numberOfDataPoints);
        }

        private boolean checkPeriodicity(int numObsPerPeriod, int numberOfDataPoints) {
            if (numObsPerPeriod == -1)
                throw new IllegalArgumentException("Must specify periodicity of seasonal");
            if (numObsPerPeriod < 4) {
                throw new IllegalArgumentException("Periodicity (numObsPerPeriod) must be >= 4");
            }
            if (numberOfDataPoints <= 2 * numObsPerPeriod) {
                throw new IllegalArgumentException(
                        "numberOfDataPoints(total length) must contain at least 2 * Periodicity (numObsPerPeriod) points");
            }
            return true;
        }
    }

    /**
     * Decompose a time series into seasonal, trend and irregular components using \code{loess}, acronym STL.
     * A new implementation of STL.  Allows for NA values, local quadratic smoothing,  post-trend smoothing, and endpoint blending.
     * The usage is very similar to that of R's built-in \code{stl()}.
     * */
    public STLResult decompose(long[] times, double[] series) {
        if (times == null || series == null || times.length != series.length)
            throw new IllegalArgumentException("times must be same length as time series");
        int n = series.length;
        int numObsPerPeriod = config.getNumObsPerPeriod();
        double[] trend = new double[n];
        double[] seasonal = new double[n];
        double[] remainder = new double[n];

        for (int i = 0; i < n; i++) {
            trend[i] = 0.0;
            seasonal[i] = 0.0;
            remainder[i] = 0.0;
        }

        if (config.getlWindow() == -1)
            config.setlWindow(STLUtility.nextOdd(config.getNumObsPerPeriod()));
        else
            config.setlWindow(STLUtility.nextOdd(config.getlWindow()));


        if (config.getSubLabels() == null) {
            int[] idx = new int[n];
            for (int i = 0; i < n; ++i)
                idx[i] = i % numObsPerPeriod + 1 ;
            config.setSubLabels(idx);
        }

        config.setsWindow(10 * n + 1);
        config.setsDegree(0);
        config.setsJump((int)Math.ceil(config.getsWindow() / 10.0));

        if (config.gettWindow() == -1) {
            /** Or use  t.window <- nextodd(ceiling(1.5 * n.p/(1 - 1.5 / s.window))) */
            config.settWindow(STLUtility.getTWindow(config.gettDegree(), config.getsDegree(), config.getsWindow(), numObsPerPeriod, config.getCritFreq()));
        }

        if (config.getsJump() == -1) config.setsJump((int)Math.ceil((double)config.getsWindow() / 10.0));
        if (config.gettJump() == -1) config.settJump((int)Math.ceil((double)config.gettWindow() / 10.0));
        if (config.getlJump() == -1) config.setlJump((int)Math.ceil((double)config.getlWindow() / 10.0));

        /** start and end indices for after adding in extra n.p before and after */
        int startIdx = numObsPerPeriod, endIdx = n - 1 + numObsPerPeriod;

        /** cycleSubIndices will keep track of what part of the
        # seasonal each observation belongs to */
        int[] cycleSubIndices = new int[n];
        double[] weight = new double[n];
        for (int i = 0; i < n; ++i) {
            cycleSubIndices[i] = i % numObsPerPeriod + 1;
            weight[i] = 1.0;
        }
        // subLabels !!
        int lenC = n + 2 * numObsPerPeriod;
        double[] C = new double[lenC];
        double[] D = new double[n];
        double[] detrend = new double[n];
        int tempSize = (int)Math.ceil((double)n / (double)numObsPerPeriod) / 2;
        ArrayList<Double> cycleSub = new ArrayList<Double>(tempSize), subWeights = new ArrayList<Double>(tempSize);
        int[] cs1 = new int[numObsPerPeriod], cs2 = new int[numObsPerPeriod];
        for (int i = 0; i < numObsPerPeriod; ++i) {
            cs1[i] = cycleSubIndices[i];
            cs2[i] = cycleSubIndices[n - numObsPerPeriod + i];
        }

        double[] ma3, L = new double[n];
        int ljump = config.getlJump(), tjump = config.gettJump();
        int lenLev = (int)Math.ceil((double)n / (double)ljump), lenTev = (int)Math.ceil((double)n / (double)tjump);
        int[] lEv = new int[lenLev], tEv = new int[lenTev];
        double weightMeanAns = 0.0;

        for (int oIter = 1; oIter <= config.getNumberOfRobustnessIterations(); ++oIter) {
            for (int iIter = 1; iIter <= config.getNumberOfInnerLoopPasses(); ++iIter) {
                /** Step 1: detrending */
                for (int i = 0; i < n; ++i)
                    detrend[i] = series[i] - trend[i];

                /** Step 2: smoothing of cycle-subseries */
                for (int i = 0; i < numObsPerPeriod; ++i) {
                    cycleSub.clear(); subWeights.clear();
                    for (int j = i; j < n; j += numObsPerPeriod) {
                        if (cycleSubIndices[j] == i + 1) {
                            cycleSub.add(detrend[j]);
                            subWeights.add(weight[j]);
                        }
                    }
                    /**
                     C[c(cs1, cycleSubIndices, cs2) == i] <- rep(weighted.mean(cycleSub,
                     w = w[cycleSubIndices == i], na.rm = TRUE), cycleSub.length + 2)
                     */
                    weightMeanAns = weightMean(cycleSub, subWeights);
                    for (int j = i; j < numObsPerPeriod; j += numObsPerPeriod)
                        if (cs1[j] == i + 1)
                            C[j] = weightMeanAns;
                    for (int j = i; j < n; j += numObsPerPeriod)
                        if (cycleSubIndices[j] == i + 1)
                            C[j + numObsPerPeriod] = weightMeanAns;
                    for (int j = 0; j < numObsPerPeriod; ++j)
                        if (cs2[j] == i + 1)
                            C[j + numObsPerPeriod + n] = weightMeanAns;
                }

                /** Step 3: Low-pass filtering of collection of all the cycle-subseries
                 # moving averages*/
                ma3 = STLUtility.cMa(C, numObsPerPeriod);

                for (int i = 0, j = 0; i < lenLev; ++i, j += ljump)
                    lEv[i] = j + 1;
                if (lEv[lenLev - 1] != n) {
                    int[] tempLev = new int[lenLev + 1];
                    System.arraycopy(lEv, 0, tempLev, 0, lenLev);
                    tempLev[lenLev] = n;
                    L = STLUtility.loessSTL(null, ma3, config.getlWindow(), config.getlDegree(), tempLev, weight, config.getlJump());
                } else {
                    L = STLUtility.loessSTL(null, ma3, config.getlWindow(), config.getlDegree(), lEv, weight, config.getlJump());
                }

                /** Step 4: Detrend smoothed cycle-subseries */
                /** Step 5: Deseasonalize */
                for (int i = 0; i < n; ++i) {
                    seasonal[i] = C[startIdx + i] - L[i];
                    D[i] = series[i] - seasonal[i];
                }

                /** Step 6: Trend Smoothing */
                for (int i = 0, j = 0; i < lenTev; ++i, j += tjump)
                    tEv[i] = j + 1;
                if (tEv[lenTev - 1] != n) {
                    int[] tempTev = new int[lenTev + 1];
                    System.arraycopy(tEv, 0, tempTev, 0, lenTev);
                    tempTev[lenTev] = n;
                    trend = STLUtility.loessSTL(null, D, config.gettWindow(), config.gettDegree(), tempTev, weight, config.gettJump());
                } else {
                    trend = STLUtility.loessSTL(null, D, config.gettWindow(), config.gettDegree(), tEv, weight, config.gettJump());
                }
            }

        }
        // Calculate remainder
        for (int i = 0; i < n; i++) {
            remainder[i] = series[i] - trend[i] - seasonal[i];
        }
        return new STLResult(trend, seasonal, remainder);
    }

    private double weightMean(ArrayList<Double> x, ArrayList<Double> w) {
        double sum = 0.0, sumW = 0.0;
        int len = x.size();
        for (int i = 0; i < len; ++i) {
            if (!Double.isNaN(x.get(i))) {
                sum += (x.get(i) * w.get(i));
                sumW += w.get(i);
            }
        }
        return sum / sumW;
    }
}
