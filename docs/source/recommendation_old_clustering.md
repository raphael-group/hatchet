# Analyze global clustering

The global clustering performed along the genome and jointly across samples is a crucial feature of HATCHet and the quality of the final results is affected by the quality of the clustering. In particular, the global clustering is performed by HATCHet's component `cluster-bins`, whose default values for the main parameters allow to deal with most of the datasets. However, noisy or special datasets need tuning of these parameters, especially because the current version of `cluster-bins` uses [scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.mixture.BayesianGaussianMixture.html), an external implementation of a GMM with a Dirichlet process prior which is not specifically designed for sequencing data. Therefore, the user needs to validate the results and improve it to obtain best-quality results, especially when considering noisy and special datasets.

The user should validate the results of the clustering, especially in noisy or suspicious cases, through the cluster figure produced by the [command CBB of plot-bins](doc_plot_bins.md) and as explained in the available [demos](https://github.com/raphael-group/hatchet#demos). More specifically, we suggest to consider the following criterion to validate the clustering: **Every pair of clusters needs to be clearly distinct in terms of RDR or BAF in at least one sample and each cluster only contains regions with similar values of RDR and BAF in all samples**.

`cluster-bins` offers two sets of parameters that the user can tune to improve the result of the clustering: (1) `cluster-bins`-specific and (2) BNPY-specific parameters.

### `cluster-bins`-specific parameters

`cluster-bins` provides two specific parameters that can be used to improve the results of the global clustering:
- `-tR` threshold defines a RDR tolerance (default value is `-tR 0.15`)
- `-tB` threshold defines a BAF tolerance (default value is `-tB 0.04`)

These tolerances are used to merge any pair of clusters which have differences in terms of RDR and BAF always lower than the given thresholds across all samples. Intuitively these two thresolds are used to draw a rectangle (with `tR` width and `tB` height in every RDR-BAF plots) around every cluster in every sample. If two clusters are in the same rectangle in every sample, they are merged. As such, the user can estimate the RDR and BAF as the values of the thresholds by using the heat figure produced by the [command BB of plot-bins](doc_plot_bins.md). In particular, higher values allow to merge more clusters, avoid overfitting, and accomodate higher noise and variances, while smaller values allow to avoid overclustering. For example, clusters which have a disc/long shape, i.e. high variance of RDR as in the [cancer WGS demo](https://github.com/raphael-group/hatchet/blob/master/examples/demo-WGS-cancer/demo-wgs-cancer.sh) or [WES demo](https://github.com/raphael-group/hatchet/blob/master/examples/demo-WES/demo-wes.sh), will require a higher value of `tR`, e.g. `-tR 0.3` or `-tR 0.4`.

### `sklearn`-specific parameters

`cluster-bins` additionally provide 3 parameters which allow to directly control the clustering computed by the scikit-learn Gaussian Mixture Model:
- `-R`: number of restarts (default is `-R 10`) which controls the number of restarts allows the method to escape local optima. While a lower number of restarts improves the running time of `cluster-bins` process, a higher number allows it to consider more solutions and choose a better clustering.
- `-K`: this value provides the maximum possible number of clusters (default is `-K 50`). Scikit-learn tends to infer close to this maximum number `K` of clusters, so if you have relatively simple data with few apparent clusters you may want to lower `K`.
- `-c`: this value (default `1`) determines the model's relative confidence in many even clusters (large `c`) vs. an uneven distribution of points into clusters, potentially with fewer than `K` components (small `c`). For experts, this is the concentration parameter for the Dirichlet process prior. Reducing this value by many orders of magnitude (e.g., `-c 0.001` or `-c 0.00001`) may enable scikit-learn to infer fewer than `K` clusters.

Note that that these clustering parameters determine the clustering *before* the merging heuristic described above is applied; if you are unsure which set of parameters to tune, consider turning off the merging heuristic (i.e., setting `-tR 0 -tB 0`) to inspect the clustering without merging.
