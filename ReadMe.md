twitter-anomalyDetection-java
=============================

A Java implementation of [twitter anomaly detection](https://github.com/twitter/AnomalyDetection).
But optimized the time complexity from o(n^2) to o(nlogn).

References
=============================
	[1. Vallis, O., Hochenbaum, J. and Kejariwal, A., (2014) “A Novel Technique for Long-Term Anomaly Detection in the Cloud”,
	 6th USENIX Workshop on Hot Topics in Cloud Computing, Philadelphia, PA.]
	(https://www.usenix.org/system/files/conference/hotcloud14/hotcloud14-vallis.pdf)
	[2. Rosner, B., (May 1983), “Percentage Points for a Generalized ESD Many-Outlier Procedure”, Technometrics, 25(2), pp. 165-172.]
	[3. STL: A Seasonal-Trend Decomposition Procedure Based on Loess](http://www.wessa.net/download/stl.pdf)

More
============================
- STL-java reference (https://github.com/brandtg/stl-java), but we implements the stl in java as [stlplus](https://github.com/hafen/stlplus) described, faster and can handle NA values (some data used stl-java will throw some exception).
- Twitter-anomalyDetection (https://github.com/twitter/AnomalyDetection), we optimize the algorithm from o(n^2) to o(nlogn)

for more information please read code, the ReadMe will update later!
