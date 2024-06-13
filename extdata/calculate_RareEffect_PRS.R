# Calculate RareEffect PRS for each individual

## Load packages and sources

library(SAIGE)
library(data.table)
library(Matrix)
library(optparse)

option_list <- list(
    make_option(c("--effectFile"), type="character", default="",
        help="Path to effect file from RareEffect main function"),
    make_option(c("--bedFile"), type="character", default="",
        help="Path to bed file containing rare variants"),
    make_option(c("--bimFile"), type="character", default="",
        help="Path to bim file containing rare variants"),
    make_option(c("--famFile"), type="character", default="",
        help="Path to fam file containing rare variants"),
    make_option(c("--groupFile"), type="character", default="",
        help="Path to group file (containing functional annotation of variants)"),
    make_option(c("--geneName"), type="character", default="",
        help="Gene name to analyze"),
    make_option(c("--variantListFile"), type="character", default="",
        help="Path to file containing list of variants to exclude (non-ultra-rare variants)"),
    make_option(c("--outputPrefix"), type="character", default="",
        help="Path to save output (without extension and gene name)")
)

## list of options
parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = FALSE)
print(args)

## Assign varibles
bedFile <- args$bedFile
bimFile <- args$bimFile
famFile <- args$famFile
effectFile <- args$effectFile
groupFile <- args$groupFile
geneName <- args$geneName

## Set PLINK object
bim <- fread(bimFile)
fam <- fread(famFile)
sampleID <- as.character(fam$V2)
modglmm <- list()
modglmm$sampleID <- sampleID
n_samples <- length(sampleID)
var_by_func_anno <- read_groupfile(groupFile, geneName)

objGeno <- SAIGE::setGenoInput(bgenFile = "",
        bgenFileIndex = "",
        vcfFile = "",
        vcfFileIndex = "",
        vcfField = "",
        savFile = "",
        savFileIndex = "",
        sampleFile = "",
        bedFile=bedFile,
        bimFile=bimFile,
        famFile=famFile,
        idstoIncludeFile = "",
        rangestoIncludeFile = "",
        chrom = "",
        AlleleOrder = "alt-first",
        sampleInModel = sampleID
    )

effect <- fread(effectFile)

# LoF
if (length(var_by_func_anno[[1]]) == 0) {
    print("LoF variant does not exist.")
}
# missense
if (length(var_by_func_anno[[2]]) == 0) {
    print("missense variant does not exist.")
}
# synonymous
if (length(var_by_func_anno[[3]]) == 0) {
    print("synonymous variant does not exist.")
}

# Remove variants not in plink file
for (i in 1:3) {
    var_by_func_anno[[i]] <- var_by_func_anno[[i]][which(var_by_func_anno[[i]] %in% objGeno$markerInfo$ID)]
}
print(str(var_by_func_anno))

# Read genotype matrix
lof_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[1]], sampleID, modglmm, macThreshold = 0)
lof_mat_collapsed <- lof_mat_collapsed_all[[1]]
lof_flipped <- lof_mat_collapsed_all[[2]]

mis_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[2]], sampleID, modglmm, macThreshold = 0)
mis_mat_collapsed <- mis_mat_collapsed_all[[1]]
mis_flipped <- mis_mat_collapsed_all[[2]]

syn_mat_collapsed_all <- collapse_matrix(objGeno, var_by_func_anno[[3]], sampleID, modglmm, macThreshold = 0)
syn_mat_collapsed <- syn_mat_collapsed_all[[1]]
syn_flipped <- syn_mat_collapsed_all[[2]]
print("Genotype matrix is loaded.")

# Remove common variants
if (args$variantListFile == "") {
    print("No variant list file is provided.")
} else {
    var_list <- fread(args$variantListFile)
    common_var <- var_list$V1
    lof_mat_collapsed <- lof_mat_collapsed[!which(colnames(lof_mat_collapsed) %in% common_var),]
    mis_mat_collapsed <- mis_mat_collapsed[!which(colnames(mis_mat_collapsed) %in% common_var),]
    syn_mat_collapsed <- syn_mat_collapsed[!which(colnames(syn_mat_collapsed) %in% common_var),]
    print("Common variants are removed.")
}

# Read effect size
effect <- as.data.frame(effect)

## LoF
if (length(which(effect$variant == "lof_UR")) == 0) {
    effect_lof_UR <- 0
} else {
    effect_lof_UR <- effect[which(effect$variant == "lof_UR"), ]$effect
}

effect_lof <- rep(effect_lof_UR, ncol(lof_mat_collapsed))

for (i in 1:(ncol(lof_mat_collapsed) - 1)) {
    if (length(which(effect$variant == colnames(lof_mat_collapsed)[i])) != 0) {
        effect_lof[i] <- effect[which(effect$variant == colnames(lof_mat_collapsed)[i]), ]$effect
    }
}

## mis
if (length(which(effect$variant == "mis_UR")) == 0) {
    effect_mis_UR <- 0
} else {
    effect_mis_UR <- effect[which(effect$variant == "mis_UR"), ]$effect
}

effect_mis <- rep(effect_mis_UR, ncol(mis_mat_collapsed))

for (i in 1:(ncol(mis_mat_collapsed) - 1)) {
    if (length(which(effect$variant == colnames(mis_mat_collapsed)[i])) != 0) {
        effect_mis[i] <- effect[which(effect$variant == colnames(mis_mat_collapsed)[i]), ]$effect
    }
}

## syn
if (length(which(effect$variant == "syn_UR")) == 0) {
    effect_syn_UR <- 0
} else {
    effect_syn_UR <- effect[which(effect$variant == "syn_UR"), ]$effect
}

effect_syn <- rep(effect_syn_UR, ncol(syn_mat_collapsed))

for (i in 1:(ncol(syn_mat_collapsed) - 1)) {
    if (length(which(effect$variant == colnames(syn_mat_collapsed)[i])) != 0) {
        effect_syn[i] <- effect[which(effect$variant == colnames(syn_mat_collapsed)[i]), ]$effect
    }
}

print("Effect size is loaded.")

## Calculate PRS

lof_prs <- lof_mat_collapsed %*% effect_lof
mis_prs <- mis_mat_collapsed %*% effect_mis
syn_prs <- syn_mat_collapsed %*% effect_syn
total_prs <- lof_prs + mis_prs + syn_prs
out <- cbind(rownames(total_prs), total_prs[,1])
out <- as.data.frame(out)
colnames(out) <- c("IID", "PRS")

## Save PRS
write.table(out, paste0(args$outputPrefix, "_", geneName, "_score.txt"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
write.table(out, "test.txt", quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
