#ifndef VCF_HPP
#define VCF_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <htslib/vcf.h>
#include <htslib/tbx.h>
#include <memory>
#include <unordered_set>

#include "bcftools_filter.hpp"
#include "bcftools_vcfnorm.hpp"

namespace VCF {

  struct site_t {
    std::string chrom;
    int32_t pos;
    std::string ref;
    std::string alt;
    bool operator==(const site_t& other_site) const {
        return (other_site.chrom == chrom && other_site.pos == pos && other_site.ref == ref && other_site.alt == alt);
    }
  };

  struct hash_site {
    std::size_t operator()(const site_t& site) const {
        return std::hash<std::string>()(site.chrom) ^ std::hash<int32_t>()(site.pos) ^ std::hash<std::string>()(site.ref) ^ std::hash<std::string>()(site.alt);
    }
  };
 
  class VcfClass{
  private:

   std::string chr;
   int32_t start = 0;
   int32_t end = 0;
   std::vector<int32_t> posSampleInVcf;
   std::vector<int32_t> posSampleInModel;
   uint32_t N0, N;
   std::vector<std::string> SampleInVcf;
 
   bool isVcfOpen; 

   std::string vcfFileName;
   size_t range_pos;
   size_t range_len;
   std::string vcfFileName_first;
   std::string vcfFileName_last;
   std::vector<htsFile*> vcf_file;
   std::vector<tbx_t*> vcf_index;
   std::vector<bcf_hdr_t*> vcf_header;
   size_t vcf_idx;
   bcf1_t* vcf_record;
   kstring_t vcf_string;
   hts_itr_t* vcf_iterator;
   std::string iterator_string;
   int has_line;
   std::vector<std::string> SampleInModel;
   bool multiple_vcf_files;
   bool iterate_over_all_lines;

   // for site iterator
   bool use_sites;
   std::unordered_set<site_t, hash_site> sites;
   std::vector<site_t> vcf_first_line;

   // for the (custom) filter
   std::string vcf_filter_string;
   std::vector<bcftools::filter_t*> vcf_filter;

   // for splitting multiallelic variants into biallelic ones
   bcf1_t* multiallelic_vcf_record;
   int multiallelic_idx;
   bcftools::args_t vcfnorm_args;

   void setPosSampleInVcf(size_t vcf_idx);
   void closeAll(size_t idx);
   void nextVcfFileName();
   void set_iterator();
   void move_forward_iterator();
   void copy_out_single_allele(int idx);
    
  public:

   std::string fmtField;

   VcfClass(std::string vcfFileName,
            std::string vcfField,
            std::string vcfFilterString,
            bool isSparseDosageInVcf,
            std::vector<std::string> SampleInModel);

   ~VcfClass();
 
 
   // setup VcfClass
   bool setVcfObj(size_t idx);
   //set up the iterator 
   void set_iterator(std::string& variantList);
   void set_iterator(std::string& chrom, int & beg_pd, int & end_pd);

   void move_forward_iterator(int i);

   bool check_iterator_end();
   void setPosSampleInVcf(std::vector<std::string> & t_SampleInModel);
   void setIsSparseDosageInVcf (bool t_isSparseDosageInVcf);
 
   // get t_dosage/genotypes of one marker

   void getOneMarker(
                                  std::string& t_ref,       // REF allele
                                  std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                  std::string& t_marker,    // marker ID extracted from genotype file
                                  uint32_t& t_pd,           // base position
                                  std::string& t_chr,       // chromosome
                                  double& t_altFreq,        // frequency of ALT allele
                                  double& t_altCounts,      // counts of ALT allele
                                  double& t_missingRate,    // missing rate
                                  double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                  bool t_isOutputIndexForMissing,               // if true, output index of missing genotype data
                                  std::vector<uint32_t>& t_indexForMissingforOneMarker,     // index of missing genotype data
                                  bool t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                  std::vector<uint32_t>& t_indexForNonZero,
                                  //if true, the marker has been read successfully
                                  bool & t_isBoolRead,
                                  arma::vec & dosages,
				  bool t_isImputation);


    void getSampleIDlist_vcfMatrix(size_t idx); 
    uint32_t getN0(){return N0;}
    uint32_t getN(){return N;}

    void closegenofile(){};

 };
 
}

#endif
