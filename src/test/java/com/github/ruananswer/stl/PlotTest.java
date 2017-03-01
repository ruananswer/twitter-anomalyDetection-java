/**
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.github.ruananswer.stl;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.ObjectMapper;
import org.jfree.data.time.Hour;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.io.File;

public class PlotTest {

  @Test
  public void testPlot() throws Exception {
    final ObjectMapper objectMapper = new ObjectMapper();
    final JsonNode tree = objectMapper.readTree(new File(this.getClass().getResource("/sample-timeseries.json").getFile()));
    final int n = tree.get("times").size();
    final long[] ts = new long[n];
    final double[] dts = new double[n];
    final double[] ys = new double[n];

    for (int i = 0; i < n; i++) {
      ts[i] = tree.get("times").get(i).asLong();
      dts[i] = ts[i];
      ys[i] = tree.get("series").get(i).asDouble();
    }

    final STLDecomposition.Config config = new STLDecomposition.Config();
    config.setNumberOfDataPoints(n);
    config.setNumObsPerPeriod(12);

    final STLDecomposition stl = new STLDecomposition(config);
    final STLResult res = stl.decompose(ts, ys);

    final File output = new File("seasonal.png");
    final File hourly = new File("stl-hourly.png");

    output.deleteOnExit();
    hourly.deleteOnExit();

    StlPlotter.plot(res, ys, dts, "New Title", Hour.class, hourly);
    StlPlotter.plot(res, ys, dts, output);
    StlPlotter.plot(res, ys, dts);

    Assert.assertTrue(output.exists());
    Assert.assertTrue(hourly.exists());

    final File exists = new File("stl-decomposition.png");
    exists.deleteOnExit();

    StlPlotter.plot(res, ys, dts, "Test Title");

    Assert.assertTrue(exists.exists());
  }
}
