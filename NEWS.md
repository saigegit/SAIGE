# SAIGE 1.5.3 (2026-11-16)

## Changed default: `--is_noadjCov` is now `FALSE` on the command line

`extdata/step2_SPAtests.R` now defaults to `--is_noadjCov=FALSE`, matching the R function
`SPAGMMATtest()` (which already defaulted to `FALSE`). With `FALSE`, covariates are regressed out of
each genotype (slower, no shortcut). `--is_noadjCov=TRUE` is the faster path and is correct again in
this release (see below). Consequence: `--is_fastTest=TRUE` without `--sparseGRMFile` now also needs
`--is_noadjCov=TRUE`; otherwise step 2 prints a message and turns `is_fastTest` off (the message is now only printed when
`--is_fastTest=TRUE` was requested). Results can
differ slightly from `--is_noadjCov=TRUE` (up to ~0.02 in -log10 p in our tests).

## Bug fix: deflated p-values with `--is_noadjCov=TRUE` when ALT is the major allele

**Affected versions:** 1.4.5 through 1.5.2, step 2, with `--is_noadjCov=TRUE` (the command-line default in these versions)
**when the single-variant test does not use a sparse GRM** (no `--sparseGRMFile` given in step 2)
**or when `--is_fastTest=TRUE`**. Step 2 runs that supply `--sparseGRMFile` with
`--is_fastTest=FALSE` (the default) do not use the affected code path. The bug is independent of
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
