#!/usr/bin/env Rscript
# Regression test: step 2 with --sparseGRMFile must use the FULL (symmetric) sparse GRM.
#
# Sparse GRMs written by createSparseGRM are symmetric MatrixMarket files that store one
# triangle only. Step 2 used to pass just the stored triangle to the C++ code, so the matrix it
# solved against was not symmetric (about half of the off-diagonal covariance was missing).
#
# This test fits a quantitative null model with the bundled sparse GRM, runs step 2, and compares
# the reported variance of the score statistic with an exact dense computation
#     var = g' P g * varRatio,   P = Sigma^-1 - Sigma^-1 X (X' Sigma^-1 X)^-1 X' Sigma^-1
# with Sigma = tau1 I + tau2 GRM (full symmetric matrix, relatednessCutoff applied).
#
# Usage (SAIGE must be installed):  Rscript tests/test_sparse_grm_symmetric.R [path/to/repo]
# Exits with status 1 on failure.
suppressMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args) >= 1) args[1] else ".")
inp <- file.path(repo, "extdata", "input"); out <- file.path(repo, "extdata", "output")
tmp <- tempfile("saige_sparse_test_"); dir.create(tmp)
sg  <- file.path(out, "sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx")
sh <- function(cmd) { rc <- system(cmd); if (rc != 0) stop("command failed: ", cmd) }

## step 1: quantitative trait, sparse GRM for null model and variance ratio
pre <- file.path(tmp, "step1")
sh(sprintf(paste("Rscript '%s/extdata/step1_fitNULLGLMM.R' --plinkFile='%s/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr_random1000'",
  "--sparseGRMFile='%s' --sparseGRMSampleIDFile='%s.sampleIDs.txt' --useSparseGRMtoFitNULL=TRUE --useSparseGRMforVarRatio=TRUE",
  "--phenoFile='%s/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt' --phenoCol=y_quantitative --covarColList=x1,x2 --qCovarColList=x2",
  "--sampleIDColinphenoFile=IID --traitType=quantitative --isCateVarianceRatio=TRUE --outputPrefix='%s' --nThreads=1 --LOCO=FALSE",
  "--IsOverwriteVarianceRatioFile=TRUE > '%s.log' 2>&1"), repo, inp, sg, sg, inp, pre, pre))

## step 2 with the sparse GRM
s2 <- file.path(tmp, "step2.txt")
sh(sprintf(paste("Rscript '%s/extdata/step2_SPAtests.R' --vcfFile='%s/genotype_100markers.vcf.gz' --vcfFileIndex='%s/genotype_100markers.vcf.gz.csi'",
  "--vcfField=GT --chrom=1 --minMAF=0 --minMAC=2 --GMMATmodelFile='%s.rda' --varianceRatioFile='%s.varianceRatio.txt'",
  "--SAIGEOutputFile='%s' --LOCO=FALSE --sparseGRMFile='%s' --sparseGRMSampleIDFile='%s.sampleIDs.txt' > '%s.log' 2>&1"),
  repo, inp, inp, pre, pre, s2, sg, sg, s2))
res <- read.table(s2, header = TRUE, sep = "\t")

## exact dense reference
load(paste0(pre, ".rda"))                       # modglmm
ids <- modglmm$sampleID; th <- modglmm$theta; X <- modglmm$X
g <- readMM(sg); sid <- read.table(paste0(sg, ".sampleIDs.txt"), colClasses = "character")[, 1]
gen <- as(as(g, "generalMatrix"), "CsparseMatrix")[match(ids, sid), match(ids, sid)]
gen@x[gen@x < 0.125] <- 0; gen <- drop0(gen)
Sigma <- th[1] * Diagonal(length(ids)) + th[2] * gen
Si <- solve(Sigma); SiX <- Si %*% X
P <- Si - SiX %*% solve(crossprod(X, SiX), t(SiX))
ratio <- as.numeric(read.table(paste0(pre, ".varianceRatio.txt"))[1, 1])

vcf <- file.path(inp, "genotype_100markers.vcf.gz")
hdr <- strsplit(grep("^#CHROM", readLines(gzfile(vcf), 200), value = TRUE)[1], "\t")[[1]]
v <- read.table(gzfile(vcf), comment.char = "#", sep = "\t", stringsAsFactors = FALSE, colClasses = "character")
colnames(v) <- hdr
dos <- function(x) c("0/0" = 0, "0/1" = 1, "1/0" = 1, "1/1" = 2)[x]
err <- sapply(res$POS, function(p) {
  gv <- as.numeric(dos(unlist(v[as.integer(v$POS) == p, ][1, ids]))); gv[is.na(gv)] <- mean(gv, na.rm = TRUE)
  ref <- as.numeric(crossprod(gv, P %*% gv)) * ratio
  abs(res$var[res$POS == p] / ref - 1) })
cat(sprintf("variants: %d   |var / exact - 1|: median %.3g, max %.3g\n", length(err), median(err), max(err)))
unlink(tmp, recursive = TRUE)
if (!(median(err) < 0.01 && max(err) < 0.05)) { cat("FAIL: step 2 variance does not match the exact dense reference\n"); quit(status = 1) }
cat("PASS\n")
