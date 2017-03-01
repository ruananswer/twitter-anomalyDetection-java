package com.github.ruananswer.stl;

import org.junit.Test;
import org.testng.Assert;

import java.io.IOException;

/**
 * Created by ruan on 16-4-14.
 */
public class STLUtilityTest {
    @Test
    public void TestLoess() throws IOException{
        // test loessSTL
        int n = 100;
        int[] x = new int[n];
        double[] y = new double[n];
        int[] m = new int[n];
        for (int i = 0; i < n; ++i) {
            x[i] = i + 1;
            m[i] = i + 1;
        }
        String[] lines = TestCommon.getResourceAsString("testLoess.csv").split("\n");
        int numData = lines.length;
        for (int i = 1; i < numData; i++) {
            String[] values = lines[i].split(",");
            String value = values[1];
            if (value.equals("NA")) {
                y[i - 1] = Double.NaN;
            } else {
                y[i - 1] = Double.valueOf(value);
            }
        }
        double[] res = STLUtility.loessSTL(x, y, 31, 0, m, null, (int)Math.ceil(31 / 10));
        Assert.assertEquals(res.length, y.length);
        return;
    }

}
