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

import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CombinedDomainXYPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.ClusteredXYBarRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.Minute;
import org.jfree.data.time.RegularTimePeriod;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;
import org.jfree.data.xy.XYDataset;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.awt.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.TimeZone;

public class StlPlotter {


  static void plot(final STLResult stlResult, double[] series, double[] times) throws IOException {
    StlPlotter.plot(stlResult, series, times, "Seasonal Decomposition");
  }

  static void plot(final STLResult stlResult, double[] series, double[] times, final File save) throws IOException {
    StlPlotter.plot(stlResult, series, times, "Seasonal Decomposition", Minute.class, save);
  }

  static void plot(final STLResult stlResult, double[] series, double[] times, final String title) throws IOException {
    StlPlotter.plot(stlResult, series, times, title, Minute.class, new File("stl-decomposition.png"));
  }

  static void plot(final STLResult stlResult, double[] series, double[] times, final String title, final Class<?> timePeriod, final File save) throws IOException {
    final ResultsPlot plot = new ResultsPlot(stlResult, series, times, title, timePeriod);
    ChartUtilities.saveChartAsPNG(save, plot.chart, 800, 600);
  }

  static void plotOnScreen(final STLResult stlResult, double[] series, double[] times, final String title) {
    final ResultsPlot plot = new ResultsPlot(stlResult, series, times, title, Minute.class);
    plot.pack();
    RefineryUtilities.centerFrameOnScreen(plot);
    plot.setVisible(true);
  }

  private static class ResultsPlot extends ApplicationFrame {

    private static final long serialVersionUID = 1L;
    private final JFreeChart chart;
    private final ChartPanel chartPanel;
    private final String title;

    private final Class<?> timePeriod;
    private final double[] series;
    private final double[] seasonal;
    private final double[] trend;
    private final double[] times;
    private final double[] remainder;

    ResultsPlot(final STLResult stlResults, double[] series, double[] times, final String title, final Class<?> timePeriod) {
      super(title);

      this.series = series;
      this.seasonal = stlResults.getSeasonal();
      this.trend = stlResults.getTrend();
      this.times = times;
      this.remainder = stlResults.getRemainder();

      this.timePeriod = timePeriod;
      this.title = title;

      this.chart = createChart();
      this.chart.removeLegend();

      this.chartPanel = new ChartPanel(chart, true, true, true, true, true);
      this.chartPanel.setPreferredSize(new Dimension(1000, 500));

      setContentPane(this.chartPanel);
    }

    private JFreeChart createChart() {

      final CombinedDomainXYPlot plot = new CombinedDomainXYPlot(new DateAxis("Time"));
      final XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer(true, false);
      final ClusteredXYBarRenderer barRenderer = new ClusteredXYBarRenderer();
      final GradientPaint black = new GradientPaint(0.0f, 0.0f, Color.black, 0.0f, 0.0f, Color.black);

      final TimeSeries seriests = new TimeSeries("Series");
      final TimeSeries seasonalts = new TimeSeries("Seasonal");
      final TimeSeries trendts = new TimeSeries("Trend");
      final TimeSeries remainderts = new TimeSeries("Remainder");

      final TimeSeries[] tsArray = new TimeSeries[]{seriests, seasonalts, trendts};
      final String[] labels = new String[]{"Series", "Seasonal", "Trend"};

      for (int i = 0; i < series.length; i++) {
        final Date d = new Date((long) times[i]);
        RegularTimePeriod rtp = RegularTimePeriod.createInstance(this.timePeriod, d, TimeZone.getDefault());
        seriests.addOrUpdate(rtp, series[i]);
        seasonalts.addOrUpdate(rtp, seasonal[i]);
        trendts.addOrUpdate(rtp, trend[i]);
        remainderts.addOrUpdate(rtp, remainder[i]);
      }

      plot.setGap(10.0);
      renderer.setSeriesPaint(0, black);
      barRenderer.setSeriesPaint(0, black);
      plot.setOrientation(PlotOrientation.VERTICAL);

      for (int i = 0; i < tsArray.length; i++) {
        final XYDataset ts = new TimeSeriesCollection(tsArray[i]);
        final XYPlot p = new XYPlot(ts, new DateAxis(labels[i]), new NumberAxis(labels[i]), renderer);
        plot.add(p);
      }

      final XYDataset rts = new TimeSeriesCollection(remainderts);
      final XYDataset sts = new TimeSeriesCollection(seriests);
      final XYDataset tts = new TimeSeriesCollection(trendts);
      final XYPlot rplot = new XYPlot(rts, new DateAxis(), new NumberAxis("Remainder"), barRenderer);
      final XYPlot seriesAndTrend = new XYPlot(sts, new DateAxis(), new NumberAxis("S & T"), renderer);

      seriesAndTrend.setDataset(1, tts);
      seriesAndTrend.setRenderer(1, renderer);

      plot.add(rplot);
      plot.add(seriesAndTrend);

      return new JFreeChart(this.title, JFreeChart.DEFAULT_TITLE_FONT, plot, true);
    }
  }

  public static void main(String[] args) throws Exception {
    List<Double> times = new ArrayList<Double>();
    List<Double> series = new ArrayList<Double>();
    List<Double> trend = new ArrayList<Double>();
    List<Double> seasonal = new ArrayList<Double>();
    List<Double> remainder = new ArrayList<Double>();

    // Read from STDIN
    String line;
    BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
    while ((line = reader.readLine()) != null) {
      String[] tokens = line.split(",");
      times.add(Double.valueOf(tokens[0]));
      series.add(Double.valueOf(tokens[1]));
      trend.add(Double.valueOf(tokens[2]));
      seasonal.add(Double.valueOf(tokens[3]));
      remainder.add(Double.valueOf(tokens[4]));
    }


    STLResult res = new STLResult(
        convert(trend),
        convert(seasonal),
        convert(remainder));

    double[] tmpSeries = convert(series);
    double[] tmpTimes = convert(times);

    if (args.length == 1) {
      plot(res, tmpSeries, tmpTimes, new File(args[0]));
    } else {
      plotOnScreen(res, tmpSeries, tmpTimes, "Seasonal Decomposition");
    }
  }

  private static double[] convert(List<Double> list) {
    double[] array = new double[list.size()];
    for (int i = 0; i < list.size(); i++) {
      array[i] = list.get(i);
    }
    return array;
  }
}
