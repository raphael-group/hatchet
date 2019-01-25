# Analyze global clustering

The global clustering performed along the genome and jointly across samples is a crucial feature of HATCHet and the quality of the final results is affected by the quality of the clustering. In particular, the global clustering is performed by HATCHet's component `cluBB`, whose default values for the main parameters allow to deal with most of the datasets. However, noisy or special datasets need tuning of these parameters, especially because the current version of `cluBB` uses [BNPY](https://bitbucket.org/michaelchughes/bnpy-dev/src/master/), an external implementation of a Dirichlet process which is not specific for sequencing data. Therefore, the user needs to validate the results and improve it to obtain best-quality results, especially when considering noisy and special datasets. 

The user should validate the results of the clustering, especially in noisy or suspicious cases, through the cluster figure produced by the [command CBB of BBot](https://github.com/raphael-group/hatchet/blob/master/doc/doc_bbot.md) and as explained in the available [demos](https://github.com/raphael-group/hatchet#demos). More specifically, we suggest to consider the following criterion to validate the clustering: **Every pair of clusters needs to be clearly distinct in terms of RDR or BAF in at least one sample and each cluster only contains regions with similar values of RDR and BAF in all samples**.

`cluBB` offers two sets of parameters that the user can tune to improve the result of the clustering: (1) `cluBB`-specific and (2) BNPY-specific parameters.

### `cluBB`-specific parameters

`cluBB` provides two specific parameters that can be used to improve the results of the global clustering:
- `-tR` threshold defines a RDR tolerance (default value is `-tR 0.15`)
- `-tB` threshold defines a BAF tolerance (default value is `-tB 0.04`)

These tolerances are used to merge any pair of clusters which have differences in terms of RDR and BAF always lower than the given thresholds across all samples. Intuitively these two thresolds are used to draw a rectangle (with `tR` width and `tB` height in every RDR-BAF plots) around every cluster in every sample. If two clusters are in the same rectangle in every sample, they are merged. As such, the user can estimate the RDR and BAF as the values of the thresholds by using the heat figure produced by the [command BB of BBot](https://github.com/raphael-group/hatchet/blob/master/doc/doc_bbot.md). In particular, higher values allow to merge more clusters, avoid overfitting, and accomodate higher noise and variances, while smaller values allow to avoid overclustering. For example, clusters which have a disc/long shape, i.e. high variance of RDR as in the [cancer WGS demo](https://github.com/raphael-group/hatchet/blob/master/examples/demo-WGS-cancer/demo-wgs-cancer.sh) or [WES demo](https://github.com/raphael-group/hatchet/blob/master/examples/demo-WES/demo-wes.sh), will require a higher value of `tR`, e.g. `-tR 0.3` or `-tR 0.4`.

### `BNPY`-specific parameters

`cluBB` additionally provide 3 parameters which allow to directly control the clustering computed by BNPY:
- `-R`: number of restarts (default is `-R 10`) which control the number of restarts for the Dirichlet process and allow to escape local optima. While a lower number of restarts improve the running time of `cluBB` process, a higher number allow to consider more solutions and choose better clustering.
- `-sf`: this value defines confidence on the size of the cluster (default is `-sf 0.01`). While lower values (e.g. `-sf 0.001` or `0.0001`) allow to consider smaller clusters, higher values (e.g. `0.1`, `1.0`) allow to consider larger clusters.
- `-K`: this value provides an initial over estimation for the number of clusters (default is `-K 15`). While lower values (e.g. `-K 10` or `-K 5`) allow to consider fewer clusters, higher values (e.g. `-K 30`, `-K 5-`) allow to consider more clusters.

As BNPY is not specific for sequencing data, it may report outlying/local/unlikely clustering solutions, especially when considering noisy or low-points datasets. In these case, the user may be able to improve the clustering by varying the values of these parameters, accordingly to the criteria previously described.
