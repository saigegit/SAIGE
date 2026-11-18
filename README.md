> **Important (2026-11-17): SAIGE 1.5.3 fixes three bugs** (details in [NEWS.md](NEWS.md)): (1) deflated step 2 p-values at ALT-major variants with `--is_noadjCov=TRUE` in 1.4.5-1.5.2 when no sparse GRM is used or with `--is_fastTest=TRUE` (#197, #206; workaround for older versions: `--is_noadjCov=FALSE`); (2) step 2 with `--sparseGRMFile` used only half of the sparse GRM; (3) an uninitialised variable in the step 1 PLINK genotype packing (#208). Please update and re-run step 2 (and step 1 for PLINK input if variance ratios matter).

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
- PLINK (bed, bim, fam), PGEN, BGEN, VCF, BCF, SAV


