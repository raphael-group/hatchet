# Improve running time

This section includes a collection of several tips for improving the overall running time.

## 1. SNP calling from known database

HATCHet allows to provide to count-alleles a list of known germline SNPs. This allows to significantly improve the performance. However, running count-alleles without this list (as by default behaviour) results in count-alleles calling germline SNPs along the whole genome and identifying more SNPs (especially including private and rare germline SNPs), which could result in a higher total number of SNPs improving quality of BAF estimations. The user can consider this trade-off.
