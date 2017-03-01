package com.github.ruananswer.stl;

/**
 * Created by on 16-4-14.
 */
public class STLUtility {
    private STLUtility() {}

    public static int nextOdd(double x) {
        int xx = (int) Math.round(x);
        return (xx%2 == 0) ? (xx + 1) : xx;
    }

    public static boolean degreeCheck(int x) {
        return (x == 0 || x == 1 || x == 2);
    }

    public static int getTWindow(int tDg, int sDg, int nS, int nP, double omega) {
        if (tDg == 0) tDg = 1;
        if (sDg == 0) sDg = 1;

        double[][] coefsA = {{0.000103350651767650, 3.81086166990428e-6}, {-0.000216653946625270, 0.000708495976681902}};
        double[][] coefsB = {{1.42686036792937, 2.24089552678906}, {-3.1503819836694, -3.30435316073732}, {5.07481807116087, 5.08099438760489}};
        double[][] coefsC = {{1.66534145060448, 2.33114333880815}, {-3.87719398039131, -1.8314816166323}, {6.46952900183769, 1.85431548427732}};
        // estimate critical frequency for seasonal
        double betac0 = coefsA[0][0] + coefsA[1][0] * omega;
        double betac1 = coefsB[0][0] + coefsB[1][0] * omega + coefsB[2][0] * omega * omega;
        double betac2 = coefsC[0][0] + coefsC[1][0] * omega + coefsC[2][0] * omega * omega;

        double fC = (1.0 - (betac0 + betac1 / nS + betac2 / (nS * nS))) / nP;

        double betat0 = coefsA[0][0] + coefsA[1][0] * omega;
        double betat1 = coefsB[0][0] + coefsB[1][0] * omega + coefsB[2][0] * omega * omega;
        double betat2 = coefsC[0][0] + coefsC[1][0] * omega + coefsC[2][0] * omega * omega;

        double betat00 = betat0 - fC;

        return nextOdd((-betat1 - Math.sqrt(betat1 * betat1 - 4.0 * betat00 * betat2)) / (2.0 * betat00));
    }

    public static double[] cMa(double[] x, int nP) {
        int i, n = x.length;
        int nn = n - nP * 2;
        int nn_p = nP;
        double ma_tmp = 0.0;

        double[] ans = new double[n - 2 * nn_p];
        double[] ma = new double[nn + nn_p + 1];
        double[] ma2 = new double[nn + 2];
        double[] ma3 = new double[nn];

        ma_tmp = 0;
        for(i = 0; i < nn_p; ++i) {
            ma_tmp = ma_tmp + x[i];
        }
        ma[0] = ma_tmp / nn_p;
        for(i = nn_p; i < nn + 2 * nn_p; ++i) {
            ma_tmp = ma_tmp - x[i - nn_p] + x[i];
            ma[i - nn_p + 1] = ma_tmp / nn_p;
        }

        ma_tmp = 0;
        for(i = 0; i < nn_p; ++i) {
            ma_tmp = ma_tmp + ma[i];
        }
        ma2[0] = ma_tmp / nn_p;

        for(i = nn_p; i < nn + nn_p + 1; ++i) {
            ma_tmp = ma_tmp - ma[i - nn_p] + ma[i];
            ma2[i - nn_p + 1] = ma_tmp / nn_p;
        }

        ma_tmp = 0;

        for(i = 0; i < 3; ++i) {
            ma_tmp = ma_tmp + ma2[i];
        }
        ans[0] = ma_tmp / 3;

        for(i = 3; i < nn + 2; ++i) {
            ma_tmp = ma_tmp - ma2[i - 3] + ma2[i];
            ans[i - 2] = ma_tmp / 3;
        }
        return ans;
    }

    public static double[] loessSTL(int[] x, double[] y, int span, int degree, int[] m, double[] weights, int jump) {
        int n = y.length;
        if (x == null) {
            x = new int[n];
            for (int i = 0; i < n; ++i)
                x[i] = i + 1;
        }
        if (weights == null) {
            weights = new double[n];
            for (int i = 0; i < n; ++i)
                weights[i] = 1.0;
        }
        int lenM = m.length;

        if (span % 2 == 0)
            span += 1;
        int s2 = (span + 1) / 2;
        int[] lIdx = new int[lenM];
        int[] rIdx = new int[lenM];
        if (n < span) {
            for (int i = 0; i < lenM; ++i) {
                lIdx[i] = 1;
                rIdx[i] = n;
            }
        } else {
            int countSmallThanS2 = 0, countLargeThanS2 = 0, countLargeThanNMinusS2 = 0;
            for (int i = 0; i < lenM; ++i) {
                if (m[i] < s2)
                    ++countSmallThanS2;
                else if (m[i] >= s2 && m[i] <= n - s2)
                    ++countLargeThanS2;
                else
                    ++countLargeThanNMinusS2;
            }
            int i = 0;
            for (; i < countSmallThanS2; ++i)
                lIdx[i] = 0;
            for (; i < countSmallThanS2 + countLargeThanS2; ++i)
                lIdx[i] = m[i] - s2;
            for (; i < countSmallThanS2 + countLargeThanS2 + countLargeThanNMinusS2; ++i)
                lIdx[i] = n - span;
            for (i = 0; i < lenM; ++i)
                rIdx[i] = lIdx[i] + span - 1;
        }
        /**
         *    aa <- abs(m - x[l_idx])
         bb <- abs(x[r_idx] - m)
         max_dist <- ifelse(aa > bb, aa, bb)
         */
        int[] maxDist = new int[lenM];
        int aa = 0, bb = 0;
        for (int i = 0; i < lenM; ++i) {
            aa = Math.abs(m[i] - x[lIdx[i]]);
            bb = Math.abs(x[rIdx[i]] - m[i]);
            maxDist[i] = (aa > bb) ? aa : bb;
        }

        if (span > n) {
            for (int i = 0; i < maxDist.length; ++i) {
                maxDist[i] = maxDist[i] + (span - n) / 2;
            }
        }
        twoArray out = cLoess(x, y, degree, span, weights, m, lIdx, maxDist);

        // do interpolation
        double[] res = new double[out.result.length];
        System.arraycopy(out.result, 0, res, 0, out.result.length);
        double[] at = new double[n];
        for (int i = 0; i < n; ++i)
            at[i] = i + 1;
        if (jump > 1) {
            res = cInterp(m, out.result, out.slopes, at);
        }
        return res;
    }

    public static class twoArray {
        public double[] result, slopes;
        public twoArray(double[] re, double[] sl) {
            result = re;
            slopes = sl;
        }
    }

    private static twoArray cLoess(int[] xx, double[] yy, int degree, int span, double[] ww, int[] m, int[] l_idx, int[] max_dist) {
        int span2, span3, offset;
        int i, j;
        double r, tmp1, tmp2;

        int n = xx.length;
        int n_m = m.length;

        double[] x = new double[span];
        double[] w = new double[span];
        double[] xw = new double[span];
        double[] x2w = new double[span];
        double[] x3w = new double[span];

        double[] result = new double[n_m];
        double[] slopes = new double[n_m];

        // variables for storing determinant intermediate values
        double a, b, c, d, e, a1, b1, c1, a2, b2, c2, det;

        span3 = span;
        if(span > n) {
            span = n;
        }

        span2 = (span - 1) / 2;

        // want to start storing results at index 0, corresponding to the lowest m
        offset = m[0];

        // loop through all values of m
        for(i = 0; i < n_m; i++) {
            a = 0.0;

            // get weights, x, and a
            for(j = 0; j < span; j++) {
                w[j] = 0.0;
                x[j] = xx[l_idx[i] + j] - (double)m[i];

                // r = std::fabs(x[j]);
                r = (x[j] > 0) ? x[j] : -x[j];
                // tricube
                tmp1 = r / max_dist[i];
                // manual multiplication is much faster than pow()
                tmp2 = 1.0 - tmp1 * tmp1 * tmp1;
                w[j] = tmp2 * tmp2 * tmp2;

                // scale by user-defined weights
                w[j] = w[j] * ww[l_idx[i] + j];

                a = a + w[j];
            }

            if(degree == 0) {
                // TODO: make sure denominator is not 0
                a1 = 1 / a;
                for(j = 0; j < span; j++) {
                    // l_i[j] = w[j] * a1;
                    result[i] = result[i] + w[j] * a1 * yy[l_idx[i] + j];
                }
            } else {
                // get xw, x2w, b, c for degree 1 or 2
                b = 0.0;
                c = 0.0;
                for(j = 0; j < span; j++) {
                    xw[j] = x[j] * w[j];
                    x2w[j] = x[j] * xw[j];
                    b = b + xw[j];
                    c = c + x2w[j];
                }
                if(degree == 1) {
                    // TODO: make sure denominator is not 0
                    det = 1 / (a * c - b * b);
                    a1 = c * det;
                    b1 = -b * det;
                    c1 = a * det;
                    for(j=0; j < span; j++) {
                        result[i] = result[i] + (w[j] * a1 + xw[j] * b1) * yy[l_idx[i] + j];
                        slopes[i] = slopes[i] + (w[j] * b1 + xw[j] * c1) * yy[l_idx[i] + j];
                    }
                } else {
                    // TODO: make sure degree > 2 cannot be specified (and < 0 for that matter)
                    // get x3w, d, and e for degree 2
                    d = 0.0;
                    e = 0.0;
                    for(j = 0; j < span; j++) {
                        x3w[j] = x[j] * x2w[j];
                        d = d + x3w[j];
                        e = e + x3w[j] * x[j];
                    }
                    a1 = e * c - d * d;
                    b1 = c * d - e * b;
                    c1 = b * d - c * c;
                    a2 = c * d - e * b;
                    b2 = e * a - c * c;
                    c2 = b * c - d * a;
                    // TODO: make sure denominator is not 0
                    det = 1 / (a * a1 + b * b1 + c * c1);
                    a1 = a1 * det;
                    b1 = b1 * det;
                    c1 = c1 * det;
                    a2 = a2 * det;
                    b2 = b2 * det;
                    c2 = c2 * det;
                    for(j=0; j < span; j++) {
                        result[i] = result[i] + (w[j] * a1 + xw[j] * b1 + x2w[j] * c1) * yy[l_idx[i] + j];
                        slopes[i] = slopes[i] + (w[j] * a2 + xw[j] * b2 + x2w[j] * c2) * yy[l_idx[i] + j];
                    }
                }
            }
        }
        return new twoArray(result, slopes);
    }

    private static double[] cInterp(int[] m, double[] fits, double[] slopes, double[] at) {
        int i, j;
        double u, h, u2, u3;

        int n_at = at.length;
        double[] ans = new double[n_at];

        j = 0; // index of leftmost vertex
        for(i = 0; i < n_at; ++i) {
            if(at[i] > m[j + 1]) j++;
            h = (m[j + 1] - m[j]);
            u = (at[i] -  m[j]) / h;
            u2 = u * u;
            u3 = u2 * u;
            ans[i] = (2 * u3 - 3 * u2 + 1) * fits[j] +
                    (3 * u2 - 2 * u3)     * fits[j + 1] +
                    (u3 - 2 * u2 + u)     * slopes[j] * h +
                    (u3 - u2)             * slopes[j + 1] * h;
        }
        return ans;
    }

}
