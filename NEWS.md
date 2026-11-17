# SAIGE 1.5.3 (2026-11-17)

This release fixes three bugs that can give wrong step 1 / step 2 results. Please read the
"who is affected" notes below.

| # | Bug | Where | Affected versions | Reported |
|---|-----|-------|-------------------|----------|
| 1 | Deflated p-values at variants where ALT is the major allele | step 2, `--is_noadjCov=TRUE` (no sparse GRM, or `--is_fastTest=TRUE`) | 1.4.5 - 1.5.2 | #197, #206 |
| 2 | Only half of the sparse GRM used | step 2 with `--sparseGRMFile` (single-variant and set-based tests) | verified on 1.5.2; the same code is in v1.3.3 and the loader in v1.0.0, so earlier versions are very likely affected (not tested) | this release |
| 3 | Uninitialised variable in the step 1 PLINK genotype packing | step 1 (full GRM and variance ratio markers), PLINK input | present since v1.0.0 (source inspection) | #208 |

## Changed default: `--is_noadjCov` is now `FALSE` on the command line

`extdata/step2_SPAtests.R` now defaults to `--is_noadjCov=FALSE`, matching the R function
`SPAGMMATtest()` (which already defaulted to `FALSE`). With `FALSE`, covariates are regressed out of
each genotype (slower, no shortcut). `--is_noadjCov=TRUE` is the faster path and is correct again in
this release (see below). Consequence: `--is_fastTest=TRUE` without `--sparseGRMFile` now also needs
`--is_noadjCov=TRUE`; otherwise step 2 prints a message and turns `is_fastTest` off (the message is now only printed when
`--is_fastTest=TRUE` was requested). Results can
differ slightly from `--is_noadjCov=TRUE` (up to ~0.02 in -log10 p in our tests).

## Bug fix 1: deflated p-values with `--is_noadjCov=TRUE` when ALT is the major allele

**Affected versions:** 1.4.5 through 1.5.2, step 2, with `--is_noadjCov=TRUE` (the command-line default in these versions)
**when the single-variant test does not use a sparse GRM** (no `--sparseGRMFile` given in step 2)
**or when `--is_fastTest=TRUE`**. Step 2 runs that supply `--sparseGRMFile` with
`--is_fastTest=FALSE` (the default) do not use *this* code path (but see bug 2 below). The bug is independent of
the genotype format (VCF/BCF, BGEN, PLINK, PGEN all share the code).
SAIGE <= 1.3.x is not affected (checked v1.3.3: the code path does not exist).

`imputeGenoAndFlip()` flips genotypes to minor-allele coding but returns the ALT allele
frequency. `scoreTestFast_noadjCov()` centred the flipped genotypes on that frequency, so the
variance of the score statistic was overestimated at every variant where ALT is the major
allele. This deflated p-values (lambda < 1). All trait types are affected (binary, quantitative
and survival were verified). For binary and survival traits SPA only re-tests variants whose
(inflated) |z| exceeds `--SPAcutoff`, so most affected variants were never corrected;
quantitative traits have no SPA step and were affected across the whole p-value range.
Score statistics (Tstat) were not affected; only their variance and hence p-values, SEs and
effect-size SEs at ALT-major variants.

Fix: `scoreTestFast_noadjCov()` now centres on the mean of the genotype vector actually tested.
Added `tests/test_allele_swap.R`, a regression test that requires identical results when REF
and ALT are swapped.

**What to do if you used 1.4.5-1.5.2 in an affected configuration:** check the step 2 log (`isSparseGRM`, `is_noadjCov`, `is_fastTest`) and, if affected, re-run step 2 with 1.5.3. Step 1 does not need to be
re-run. Workaround on older versions: `--is_noadjCov=FALSE`.

Reported in #197 and #206; thanks to Karen Oliver and Melanie Bahlo (Walter and Eliza Hall
Institute / University of Melbourne) for the detailed analysis.

## Bug fix 2: step 2 with `--sparseGRMFile` used only half of the sparse GRM

Sparse GRMs written by `createSparseGRM` are symmetric MatrixMarket files, which store only one
triangle. `setSparseSigma_new()` passed just those stored entries to the C++ code, which builds
Sigma from them as-is and solves against it. The matrix used was therefore not symmetric: about half
of the off-diagonal (relatedness) covariance was missing, which changes the variance of the score
statistic (and so p-values and SEs) in single-variant and set-based tests that use the sparse GRM.

Fix: expand to general storage (both triangles) after reading the file. New regression test
`tests/test_sparse_grm_symmetric.R` compares the step 2 variance with an exact dense computation.
On the bundled example data (quantitative trait, sparse GRM) the median relative error of the
variance was 36% (max about 400%, including a negative variance) before the fix and 0.1% (max 1%)
after. On the bundled binary example the effect was small (about 1% in the variance), so the
impact depends on the trait and on how related the cohort is.

**Who should re-run:** anyone who ran step 2 with `--sparseGRMFile` (including all set-based /
SAIGE-GENE+ runs that use a sparse GRM). Step 1 does not need to be re-run for this bug. Step 1's own
handling of the sparse GRM (`--useSparseGRMtoFitNULL`) was not examined for this release.

## Bug fix 3: uninitialised variable in the step 1 PLINK genotype packing

In `src/SAIGE_fitGLMM_fast.cpp`, the byte that packs four 2-bit genotypes (`geno2`) was declared
without an initial value, while `setGenotype()` only ORs bits in. The first packed byte of every
marker (samples 1-4) therefore started from stale stack memory; a non-carrier (0b00) could look like
a carrier. This affects the markers kept for the GRM and those kept for variance ratio estimation,
so the effect is on the null model and the variance ratios (reported effect about 5% on the
rare-variant variance ratio in #208; the size depends on the build and the data).

Fix: initialise `geno2` to 0. Verified on a PLINK set with 40 copies of one MAC-12 marker: with the
fix the variance ratios are identical with and without a compiler flag that fills uninitialised
locals with junk; without the fix they differ from the fixed result and between builds.

**Who should re-run:** results from step 1 with PLINK input may change slightly; re-fitting the null
model with 1.5.3 is recommended for studies where the variance ratios matter (for example rare
variants, categorical variance ratios).

