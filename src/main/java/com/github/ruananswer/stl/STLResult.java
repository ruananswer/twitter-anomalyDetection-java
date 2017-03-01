package com.github.ruananswer.stl;

/**
 * The STL decomposition of a time series.
 * <p>
 * getData() == getTrend() + getSeasonal() + getRemainder()
 * </p>
 */
public class STLResult {
    private final double[] trend;
    private final double[] seasonal;
    private final double[] remainder;

    public STLResult(double[] trend, double[] seasonal,
                     double[] remainder) {
        this.trend = trend;
        this.seasonal = seasonal;
        this.remainder = remainder;
    }

    public double[] getTrend() {
        return trend;
    }

    public double[] getSeasonal() {
        return seasonal;
    }

    public double[] getRemainder() {
        return remainder;
    }
}