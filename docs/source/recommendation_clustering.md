# Analyze global clustering

The global clustering performed along the genome and jointly across samples is a crucial feature of HATCHet and the quality of the final results is strongly affected by the quality of the clustering. This global clustering is performed by HATCHet's component `cluster-bins`, whose default values are suitable for many datasets. However, for ideal results on specific datasets these parameters may need to be modified.

The module `cluster-bins` incorporates genomic position to improve clustering using a Gaussian hidden Markov model (GHMM), as opposed to the position-agnostic Gaussian mixture model (GMM) used in `cluster-bins-gmm` and described in the original HATCHet publication. This page describes how to tune the parameters of `cluster-bins` -- for recommendations on `cluster-bins-gmm`, see [this page](recommendation_old_clustering.md) instead.

The user should validate the results of the clustering, especially in noisy or suspicious cases, through the cluster figures produced by [plot-bins](doc_plot_bins.html) and [plot-bins-1d2d](doc_plot_bins_1d2d.html). More specifically, we suggest the following criteria to evaluate the clustering:

1. Every pair of clusters should be clearly distinct in terms of RDR or BAF in at least one sample, and
2. Each cluster should contain regions with similar values of RDR and BAF in all samples


`cluster-bins` offers several parameters that can be used to tune the clustering.
## Number of clusters
By default, `cluster-bins` tries several possible values for the number `K` of clusters and selects the one that maximizes the silhouette score. In practice, this tends to *underestimate* the number of clusters that are visually apparent. This can be modified by

1. Setting the parameters `--minK` and `--maxK` which specify the minimum and maximum number of clusters to consider, or
2. Setting the parameter `--exactK` to fix the number of clusters to a given value.

## Model parameters
Some parameters of the model can be tuned to change the clustering. The most useful one is the value `--tau`, which corresponds to the probability of transitioning between different copy-number states (i.e., the initial value for off-diagonal entries in the transition matrix). In practice, `tau` controls the balance between global information (RDR and BAf across samples) and local information (assigning adjacent bins to the same cluster): smaller values of `tau` put more weight on *local* information, and larger values of `tau` put more weight on *global* information.
 It may be appropriate to reduce `tau` by several orders of magnitude for noisier or lower-coverage datasets, as in this case the global RDR/BAF values are less reliable.

 Other parameters (below) are available to change the structure of the model, although in practice I have not found them particularly helpful in tuning the clustering.

| Name | Description | Usage | Default |
|------|-------------|-------|---------|
| `-t`, `--transmat` | Type of transition matrix to infer | `fixed` (to off-diagonal = tau), `diag` (all diagonal elements are equal, all off-diagonal elements are equal) or `full` (freely varying) | `diag` |
| `-c`, `--covar` | Type of covariance matrix to infer | options described in [hmmlearn documentation](https://hmmlearn.readthedocs.io/en/latest/api.html#hmmlearn.hmm.GaussianHMM) | `diag` |
| `-x`, `--decoding` | Decoding algorithm to use to infer final estimates of states | `map` for MAP inference, `viterbi` for Viterbi algorithm | `map` |
