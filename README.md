HERE we are maintaining an newly improved stable version of SAIGE and SAIGE-GENE+. 
Please find the https://saigegit.github.io/SAIGE-doc/ for documentation.


SAIGE is an R package developed with Rcpp for genome-wide association tests in large-scale data sets and biobanks. The method

- accounts for sample relatedness based on the generalized mixed models
- allows for model fitting with either full or sparse genetic relationship matrix (GRM)
- works for quantitative and binary traits
- handles case-control imbalance of binary traits
- computationally efficient for large data sets
- performs single-variant association tests
- provides effect size estimation through Firth's Bias-Reduced Logistic Regression
- performs conditional association analysis

SAIGE-GENE (now known as SAIGE-GENE+) are new method extension in the R package for testing rare variant in set-based tests.
- performs BURDEN, SKAT, and SKAT-O tests
- allows for tests on multiple minor allele frequencies cutoffs and functional annotations
- allows for specifying weights for markers in the set-based tests
- performs conditional analysis to identify associations independent from nearly GWAS signals


The package takes genotype file input in the following formats
- PLINK (bed, bim, fam), BGEN, VCF, BCF, SAV

To build SAIGE using Docker - clone the repo, revert to desired version (x.x.x) if needed & run the following command from the root directory (prune to remove intermediate containers from build) :
- ```docker build -t saige:x.x.x --build-arg BASE_IMAGE=ubuntu:xx.xx . && docker image prune```
- The `--build-arg` flag is optional : `BASE_IMAGE` will default to ubuntu:24.04 LTS. For short term builds, specify the latest intermediate release for updated packages.

