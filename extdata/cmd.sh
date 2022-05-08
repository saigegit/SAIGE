#When a sparse GRM is used to fit the null model (--useSparseGRMtoFitNULL=TRUE), variance ratios are estiamted using markers randomly selected from (--plinkFile). Categorical variance ratios are estimated (--isCateVarianceRatio=TRUE)

Rscript step1_fitNULLGLMM.R     \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr  \
    --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
    --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
    --useSparseGRMtoFitNULL=TRUE    \
    --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
    --phenoCol=y_binary \
    --covarColList=x1,x2,a9,a10 \
    --qCovarColList=a9,a10  \
    --sampleIDColinphenoFile=IID \
    --traitType=binary        \
    --IsOverwriteVarianceRatioFile=TRUE	\
    --isCateVarianceRatio=TRUE      \
    --outputPrefix=./output/example_binary_sparseGRM

#
 Rscript step2_SPAtests.R        \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/genotype_100markers_bgen_groupTest_out_sparseGRMforStep1.txt \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile=./input/samplelist.txt \
     --GMMATmodelFile=./output/example_binary_sparseGRM.rda \
     --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
     --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous        \
     --maxMAF_in_groupTest=0.0001,0.001,0.01	\
     --LOCO=FALSE	\
     --varianceRatioFile=./output/example_binary_sparseGRM.varianceRatio.txt	\
     --is_fastTest=TRUE




     Rscript step2_SPAtests.R        \
        --vcfFile=./input/genotype_100markers_missGT.withchr.vcf.gz    \
        --vcfFileIndex=./input/genotype_100markers_missGT.withchr.vcf.gz.csi     \
        --vcfField=GT   \
        --SAIGEOutputFile=./output/genotype_100markers_vcf_groupTest_out.txt \
        --LOCO=FALSE    \
        --minMAF=0 \
        --minMAC=0.5 \
        --sampleFile=./input/samplelist.txt \
        --GMMATmodelFile=./output/example_binary_fullGRM.rda \
        --varianceRatioFile=./output/example_binary_cate.varianceRatio.txt      \
        --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
        --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
        --groupFile=./input/group_new_chrposa1a2_withchr.txt    \
        --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous        \
        --maxMAF_in_groupTest=0.0001,0.001,0.01	\
	--chrom=chr1
