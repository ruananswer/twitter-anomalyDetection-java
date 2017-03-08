package com.github.ruananswer.statistics;

import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.Arrays;

/**
 * Created by ruan on 16-4-22.
 */
public class QuickMedianTest {
    @Test
    public void simpleArray() {
        double[] testArray = {
                2, 1, 6, 5, 4, 3
        };
        testCorrectness(testArray);
    }

    @Test
    public void bigArray() {
        double[] testArray = new double[1000];
        for (int i = 0; i < testArray.length; i++) {
            testArray[i] = Math.random();
        }
        testCorrectness(testArray);
    }

    private void testCorrectness(double[] arr) {
        QuickMedians stat = new QuickMedians(arr);

        Arrays.sort(arr);
        double median = 0;
        if (arr.length % 2 == 0)
            median = (arr[arr.length / 2 - 1] + arr[arr.length / 2]) / 2.0;
        else
            median = arr[arr.length / 2];
        Assert.assertEquals(stat.getMedian(), median, 0.001);
    }
}
