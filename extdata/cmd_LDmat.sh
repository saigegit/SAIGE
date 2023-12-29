
#generate LD for all markers in bgen for each set in group file --annotation_in_groupTest=ALL
Rscript step3_LDmat.R	\
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/LDmat \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --sampleFile=./input/samplelist.txt \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest=ALL        \
     --maxMAF_in_groupTest=0.5


##generate LD for lof and missense markers  in bgen for each set in group file --annotation_in_groupTest="lof;missense"
Rscript step3_LDmat.R   \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/LDmat \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --sampleFile=./input/samplelist.txt \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest="lof;missense"        \
     --maxMAF_in_groupTest=0.5


###generate LD for all markers in VCF for each set in group file --annotation_in_groupTest=ALL
Rscript step3_LDmat.R   \
      --vcfFile=./input/genotype_100markers.vcf.gz    \
      --vcfFileIndex=./input/genotype_100markers.vcf.gz.csi     \
      --vcfField=GT   \
      --SAIGEOutputFile=./output/LDmat \
      --chrom=1 \
      --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest=ALL        \
     --maxMAF_in_groupTest=0.5


#generate LD for all markers in bgen for each set in group file --annotation_in_groupTest=ALL
#use a subset of samples specified in --sample_include_inLDMat_File=./input/sample.subset.txt (one column with sample IDs, no header)
/bin/time -o LDmat.runinfo.txt -v Rscript step3_LDmat.R   \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./output/LDmat \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --sampleFile=./input/samplelist.txt \
     --groupFile=./input/group_new_chrposa1a2.txt    \
     --annotation_in_groupTest=ALL        \
     --maxMAF_in_groupTest=0.5	\
     --sample_include_inLDMat_File=./input/sample.subset.txt
