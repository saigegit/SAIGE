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
    --outputPrefix=./output/example_binary_sparseGRM_temo




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




     Rscript step1_fitNULLGLMM.R     \
    --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
    --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
    --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly_22chr \
    --useSparseGRMtoFitNULL=TRUE    \
    --phenoFile=./input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
    --phenoCol=y_binary \
    --covarColList=x1,x2 \
    --qCovarColList=x2  \
    --sampleIDColinphenoFile=IID \
    --cateVarRatioMinMACVecExclude="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude="1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --traitType=binary        \
    --isCateVarianceRatio=TRUE	\
    --outputPrefix=./output/example_binary_sparseGRM	\
    --IsOverwriteVarianceRatioFile=TRUE



      Rscript step2_SPAtests.R        \	      
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/genotype_100markers_bgen_groupTest_out.txt \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile=./input/samplelist.txt \
     --GMMATmodelFile=./output/example_binary_sparseGRM.rda \
     --varianceRatioFile=output/example_binary_sparseGRM.varianceRatio.col1.txt	\
     --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
     --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt  \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest="intergenic"        \
     --maxMAF_in_groupTest=0.5 \
     --is_output_markerList_in_groupTest=TRUE \
 --LOCO=FALSE \
     --is_fastTest=FALSE	\
         --cateVarRatioMinMACVecExclude="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude="1.5,2.5,3.5,4.5,5.5,10.5,20.5"



       Rscript step2_SPAtests.R        \
    --vcfFile=./input/genotype_100markers_missGT.withchr.vcf.gz    \
    --vcfFileIndex=./input/genotype_100markers_missGT.withchr.vcf.gz.csi     \
    --vcfField=GT   \
    --SAIGEOutputFile=./output/genotype_100markers_vcf_groupTest_out.txt \
    --LOCO=FALSE    \
    --minMAF=0 \
    --minMAC=0.5 \
    --GMMATmodelFile=./output/example_binary_fullGRM.rda \
    --varianceRatioFile=output/example_binary_sparseGRM.varianceRatio.col1.txt	\
    --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
    --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
    --groupFile=./input/group_new_chrposa1a2_withchr.txt    \
    --annotation_in_groupTest=intergenic        \
    --maxMAF_in_groupTest=0.0001	\
             --cateVarRatioMinMACVecExclude="0.5,1.5,2.5,3.5,4.5,5.5,10.5,20.5" \
    --cateVarRatioMaxMACVecInclude="1.5,2.5,3.5,4.5,5.5,10.5,20.5"



Rscript step2_SPAtests.R        \
      --vcfFile=./input/genotype_100markers.vcf.gz    \
      --vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
      --vcfField=GT   \
      --SAIGEOutputFile=./output/genotype_100markers_marker_vcf_cond.txt_quantitative \
      --chrom=1       \
      --minMAF=0 \
      --minMAC=20 \
      --GMMATmodelFile=./output/example_quantitative.rda \
      --varianceRatioFile=./output/example_quantitative.varianceRatio.txt   \
      --is_output_moreDetails=TRUE    \
      --condition=1:13:A:C,1:79:A:C


 Rscript step2_SPAtests.R        \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/genotype_100markers_bgen_groupTest_out_cond.txt_quantitative \
     --chrom=1 \
     --LOCO=TRUE    \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile=./input/samplelist.txt \
     --GMMATmodelFile=./output/example_quantitative_sparseGRM.rda \
     --varianceRatioFile=./output/example_quantitative_sparseGRM.varianceRatio.txt      \
     --sparseGRMFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx   \
     --sparseGRMSampleIDFile=output/sparseGRM_relatednessCutoff_0.125_1000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt     \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous        \
     --maxMAF_in_groupTest=0.0001,0.001,0.01	\
     --condition=1:30:A:C,1:79:A:C


Rscript step2_SPAtests.R        \
      --vcfFile=./input/genotype_100markers.vcf.gz    \
      --vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
      --vcfField=GT   \
      --SAIGEOutputFile=./output/genotype_100markers_marker_vcf_binary_ER_test_temp \
      --chrom=1       \
      --minMAF=0 \
      --minMAC=1 \
      --GMMATmodelFile=./output/example_binary.rda \
      --varianceRatioFile=./output/example_binary.varianceRatio.txt   \
      --is_output_moreDetails=TRUE	\
      --rangestoIncludeFile=./input/includerange_oneline.txt


##time-to-event phenotypes
Rscript step1_fitNULLGLMM.R             --plinkFile=./input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly         --phenoFile=./input/pheno_1000samples_survival.txt         --phenoCol=casecontrol         --covarColList=X         --eventTimeCol=AgeOfEventFinal         --sampleIDColinphenoFile=IND_ID         --traitType=survival                --outputPrefix=./output/example_survival         --nThreads=4            --LOCO=FALSE            --minMAFforGRM=0.01             --skipModelFitting=FALSE    --tol=0.01  --isCovariateOffset=FALSE --IsOverwriteVarianceRatioFile=TRUE  > test.log.survival
