#!/usr/bin/env Rscript
# Regression test: single-variant step 2 results must not depend on which allele is ALT.
#
# Swaps REF/ALT (and recodes genotypes) in the bundled example VCF, runs step 1 once and
# step 2 on both versions of the VCF, and requires identical p-values / score statistics /
# variances. Guards against the --is_noadjCov bug (SAIGE 1.4.5-1.5.2) where the variance
# was inflated whenever ALT was the major allele (see GitHub issues #197, #206).
#
# Usage (SAIGE must be installed; needs bgzip and tabix on PATH):
#   Rscript tests/test_allele_swap.R [path/to/repo]
# Exits with status 1 on failure.

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(if (length(args) >= 1) args[1] else ".")
inp  <- file.path(repo, "extdata", "input")
tmp  <- tempfile("saige_swap_test_"); dir.create(tmp)
tol  <- 1e-6

for (tool in c("bgzip", "tabix"))
  if (Sys.which(tool) == "") stop("required tool not found on PATH: ", tool)

sh <- function(cmd) { rc <- system(cmd); if (rc != 0) stop("command failed: ", cmd) }

## 1. original VCF and an ALT/REF-swapped copy
orig <- file.path(tmp, "orig.vcf.gz"); swp <- file.path(tmp, "swapped.vcf.gz")
file.copy(file.path(inp, "genotype_100markers.vcf.gz"), orig)
awk <- paste0("awk 'BEGIN{FS=OFS=\"\\t\"} /^#/{print;next} {t=$4;$4=$5;$5=t; ",
              "for(i=10;i<=NF;i++){g=$i; if(g==\"0/0\")$i=\"1/1\"; else if(g==\"1/1\")$i=\"0/0\"; ",
              "else if(g==\"0|0\")$i=\"1|1\"; else if(g==\"1|1\")$i=\"0|0\"} print}'")
sh(sprintf("gzip -dc '%s' | %s | bgzip -c > '%s'", orig, awk, swp))
for (f in c(orig, swp)) sh(sprintf("tabix -f -C -p vcf '%s'", f))

## 2. step 1 (full GRM, binary trait)
pre <- file.path(tmp, "step1")
sh(sprintf(paste("Rscript '%s/extdata/step1_fitNULLGLMM.R'",
  "--plinkFile='%s/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr_random1000'",
  "--phenoFile='%s/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt'",
  "--phenoCol=y_binary --covarColList=x1,x2 --qCovarColList=x2 --sampleIDColinphenoFile=IID",
  "--traitType=binary --outputPrefix='%s' --nThreads=1 --LOCO=FALSE",
  "--IsOverwriteVarianceRatioFile=TRUE > '%s.log' 2>&1"),
  repo, inp, inp, pre, pre))

## 3. step 2 on both codings, default --is_noadjCov=TRUE
step2 <- function(vcf, out)
  sh(sprintf(paste("Rscript '%s/extdata/step2_SPAtests.R' --vcfFile='%s' --vcfFileIndex='%s.csi'",
    "--vcfField=GT --chrom=1 --minMAF=0 --minMAC=2 --GMMATmodelFile='%s.rda'",
    "--varianceRatioFile='%s.varianceRatio.txt' --SAIGEOutputFile='%s' --LOCO=FALSE",
    "--is_noadjCov=TRUE > '%s.log' 2>&1"), repo, vcf, vcf, pre, pre, out, out))
o_out <- file.path(tmp, "orig.txt"); s_out <- file.path(tmp, "swapped.txt")
step2(orig, o_out); step2(swp, s_out)

## 4. compare on shared variants
o <- read.table(o_out, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
s <- read.table(s_out, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
common <- intersect(o$POS, s$POS)
if (length(common) < 50) stop("too few shared variants: ", length(common))
o <- o[match(common, o$POS), ]; s <- s[match(common, s$POS), ]
d_logp <- max(abs(-log10(o$p.value) - -log10(s$p.value)))
d_var  <- max(abs(o$var / s$var - 1))
d_T    <- max(abs(abs(o$Tstat) - abs(s$Tstat)))
cat(sprintf("shared variants: %d\nmax |d(-log10 p)|: %.3g\nmax |var ratio - 1|: %.3g\nmax ||Tstat| diff|: %.3g\n",
            length(common), d_logp, d_var, d_T))
ok <- d_logp < tol && d_var < tol && d_T < tol
unlink(tmp, recursive = TRUE)
if (!ok) { cat("FAIL: results depend on allele coding\n"); quit(status = 1) }
cat("PASS\n")
