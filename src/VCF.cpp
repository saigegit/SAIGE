// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
 
#include "savvy/reader.hpp"
#include "savvy/region.hpp"
#include "variant_group_iterator.hpp"


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <Rcpp.h>
#include <stdlib.h>
#include <cstring>
#include <limits>
#include <map>
#include <sstream>


#include "VCF.hpp"

// Forward declarations for helper functions
void push_filter(std::vector<std::tuple<std::string, float, float>> & filters_parsed, std::tuple<std::string, float, float> & cur_parsed);
std::tuple<std::string, float, float> parse_filter(std::string & filter);
std::vector<std::tuple<std::string, float, float>> parse_filters(std::string & filters);
std::vector<std::tuple<std::string, float, float>> consolidate_filters(
  const std::vector<std::tuple<std::string, float, float>>& filters);


using namespace std;
 
namespace VCF {
 
  VcfClass::VcfClass(std::string t_vcfFileName,
            std::string t_vcfFileIndex,
            std::string t_vcfField,
            std::string t_vcfFilters,
            bool t_isSparseDosageInVcf,
            std::vector<std::string> t_SampleInModel)
   {
     bool isVcfOpen;
     if(t_vcfFilters != ""){
      m_vcfFilters = parse_filters(t_vcfFilters);
     }
     isVcfOpen = setVcfObj(t_vcfFileName,
              t_vcfFileIndex,
              t_vcfField);
 
     setPosSampleInVcf(t_SampleInModel);
     setIsSparseDosageInVcf(t_isSparseDosageInVcf); 
   }
 
  
   bool VcfClass::setVcfObj(const std::string & t_vcfFileName,
                  const std::string & t_vcfFileIndex,
                  const std::string & t_vcfField)
   {
     //using namespace Rcpp;
     m_marker_file = savvy::reader(t_vcfFileName);
     bool isVcfOpen = m_marker_file.good();
    std::string fmtField;
    if(isVcfOpen){
    fmtField = t_vcfField;
    std::cout << "Open VCF done" << std::endl;
    cout << "To read the field " << t_vcfField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << m_marker_file.headers().size() << endl;
    m_N0 = m_marker_file.samples().size();
    std::cout << "Number of samples in the vcf file: " << m_N0 << endl;

    if (std::find_if(m_marker_file.format_headers().begin(), m_marker_file.format_headers().end(),
      [&](const savvy::header_value_details& h) { return h.id == fmtField; }) == m_marker_file.format_headers().end()) {
      // Regarding dosage fields:
      // DS and GT are equivalent.
      // That is, a given VCF may list its dosage field as "GT" rather than "DS".
      if ((fmtField == "DS" || fmtField == "GT") && std::find_if(m_marker_file.format_headers().begin(), m_marker_file.format_headers().end(),
        [&](const savvy::header_value_details& h) { return h.id == "HDS"; }) != m_marker_file.format_headers().end()) {
        fmtField = "HDS";
      } else {
        std::cerr << "ERROR: vcfField (" << fmtField << ") not present in genotype file." << std::endl;
        return false;
      }
    }
    m_fmtField = fmtField;
  } else {
    std::cerr << "WARNING: Open VCF failed" << std::endl;
  }
  return(isVcfOpen);
 }
 
  void VcfClass::set_iterator(std::string & chrom, int & beg_pd, int& end_pd)
  {
      m_it_ = variant_group_iterator(m_marker_file, savvy::region(chrom, beg_pd, end_pd));
      //variant_group_iterator end{};
  }


  void VcfClass::set_iterator(std::string & variantList)
  {
      m_it_ = variant_group_iterator(m_marker_file, variantList);
      //variant_group_iterator end{}; 
  }

  void VcfClass::move_forward_iterator(int i)
  {	

	for(int j = 0; j < i; j++){
		++m_it_;
	}
  }

  bool VcfClass::check_iterator_end()
  {	 
	bool isEndFile = false;
        variant_group_iterator end{};	
        if(m_it_ != end){
		isEndFile = false;
	}else{
		isEndFile = true;
	}
	return(isEndFile);
	//std::cout << "isEndFile " << isEndFile << std::endl;	
  }

  void VcfClass::setPosSampleInVcf(std::vector<std::string> & t_SampleInModel)
  {
     std::cout << "Setting position of samples in VCF files...." << std::endl;
     m_N = t_SampleInModel.size();
   std::cout << "m_N " << m_N << std::endl;

     getSampleIDlist_vcfMatrix();
     Rcpp::CharacterVector SampleInVcf(m_N0);
     for(uint32_t i = 0; i < m_N0; i++)
       SampleInVcf(i) = m_SampleInVcf.at(i);

   Rcpp::CharacterVector SampleInModel = Rcpp::CharacterVector::create();

   if(m_N != 0){
     //Rcpp::CharacterVector SampleInModel(m_N);
     for(uint32_t i = 0; i < m_N; i++)
       SampleInModel.push_back(t_SampleInModel.at(i));	     
       //SampleInModel(i) = t_SampleInModel.at(i);
   }else{
     m_N = m_N0;
     //Rcpp::CharacterVector SampleInModel(m_N);
     for(uint32_t i = 0; i < m_N; i++)
        SampleInModel.push_back(SampleInVcf(i));
   }	   

     Rcpp::IntegerVector posSampleInVcf = Rcpp::match(SampleInModel, SampleInVcf);
     for(uint32_t i = 0; i < m_N; i++){
       if(Rcpp::IntegerVector::is_na(posSampleInVcf.at(i)))
          Rcpp::stop("At least one subject requested is not in VCF file.");
     }

     Rcpp::IntegerVector posSampleInModel = Rcpp::match(SampleInVcf, SampleInModel);
     m_posSampleInModel.resize(m_N0);
     for(uint32_t i = 0; i < m_N0; i++){
       if(Rcpp::IntegerVector::is_na(posSampleInModel.at(i))){
         m_posSampleInModel.at(i) = -1;
       }else{
         m_posSampleInModel.at(i) = posSampleInModel.at(i) - 1;   // convert "starting from 1" to "starting from 0"
       }
     }
  }

   // get dosages/genotypes of one marker   
   void VcfClass::getOneMarker(
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
			  bool t_isImputation)
			  //std::vector<double> & t_dosage) 
   {
     bool isReadVariant = true;
     variant_group_iterator end{};
     if(m_it_ != end){
       if (m_it_->alts().size() != 1)
       {
         std::cerr << "Warning: skipping multiallelic variant" << std::endl;
       }
       t_chr = m_it_->chromosome();
       //t_pd = std::to_string(m_it_->position());
       t_pd = m_it_->position();
       t_ref = m_it_->ref();
       t_alt = m_it_->alts()[0];
       t_marker = m_it_->id(); 
       float markerInfo = 1.f;
       if(t_isImputation){
       	m_it_->get_info("R2", markerInfo);
       	t_imputeInfo = double(markerInfo);
       }else{
	t_imputeInfo = 1.0;
       }

       double dosage;
       t_altCounts = 0;
       int missing_cnt = 0;
       int filter_cnt = 0;

       dosages.clear();
       dosages.set_size(m_N);
       dosages.fill(arma::fill::zeros);
       //t_dosage.clear();
       //t_dosage.reserve(m_N);

       // compressed_vector: savvy construct, vector of zeros (sparse vector)
       savvy::compressed_vector<float> variant_dosages;

       // std::vector<std::tuple<std::string, float, float>> filters; // format field name, lower bound, upper bound
       std::vector<savvy::compressed_vector<float>> format_data(m_vcfFilters.size());

       for (size_t i=0; i<m_vcfFilters.size(); i++){
          m_it_->get_format(std::get<0>(m_vcfFilters[i]), format_data[i]);
       }

       m_it_->get_format(m_fmtField, variant_dosages);
       std::size_t ploidy = m_N0 ? variant_dosages.size() / m_N0 : 1;
       savvy::stride_reduce(variant_dosages, ploidy);
       for (auto dose_it = variant_dosages.begin(); dose_it != variant_dosages.end(); ++dose_it) {
        // j: offset, brings us to next dose value
        // this "offset" corresponds to column (e.g. individual)
        int j = dose_it.offset();
        if(m_posSampleInModel[j] >= 0) {

          bool passes_filters = true;
          for (size_t i=0; i<m_vcfFilters.size(); i++){
            float value = format_data[i][j];
            float lower_bound = std::get<1>(m_vcfFilters[i]);
            float upper_bound = std::get<2>(m_vcfFilters[i]);
            if (value < lower_bound || value > upper_bound){
              passes_filters = false;
              ++filter_cnt; // book keeping as of now; not used anywhere yet
              break;
            }
          }
          if (!passes_filters || std::isnan(*dose_it)) {
            dosages[m_posSampleInModel[j]] = -1;		  
            ++missing_cnt;
            t_indexForMissingforOneMarker.push_back(m_posSampleInModel[j]);
          }else{
              dosages[m_posSampleInModel[j]] = *dose_it; 		 	 
              if(*dose_it > 0){
                      t_indexForNonZero.push_back(m_posSampleInModel[j]);
              }		    
              //dosagesforOneMarker[genetest_sample_idx_vcfDosage[i]] = *dose_it;
              t_altCounts = t_altCounts + *dose_it;
          }
        }
      }

       if(missing_cnt > 0){
         if(missing_cnt == m_N){
           t_altFreq = 0;
         }else{
           t_altFreq = t_altCounts / 2 / (double)(m_N - missing_cnt);
         }
	 t_missingRate = missing_cnt/double(m_N);
        }else{
          t_altFreq = t_altCounts / 2 / (double)(m_N);
	 t_missingRate = 0; 
        }
     }else{
       isReadVariant = false;	     
       std::cout << "Reach the end of the vcf file" << std::endl;
     }
     t_isBoolRead = isReadVariant;  
     //return(isReadVariant);    	  
   }
 
void VcfClass::setIsSparseDosageInVcf (bool t_isSparseDosageInVcf){
  m_isSparseDosagesInVcf = t_isSparseDosageInVcf;
}

void VcfClass::getSampleIDlist_vcfMatrix(){
  std::vector< std::string > sampleIDList (m_marker_file.samples().begin(), m_marker_file.samples().end());
  m_SampleInVcf = sampleIDList;
}

}

std::vector<std::tuple<std::string, float, float>> parse_filters(std::string & filters){
  std::vector<std::tuple<std::string, float, float>> filters_parsed;

  std::stringstream ss(filters);
  std::string filter;
  
  while(std::getline(ss, filter, ',')){
    auto cur_parsed = parse_filter(filter);
    push_filter(filters_parsed, cur_parsed);
  }
  return filters_parsed;
}

void push_filter(std::vector<std::tuple<std::string, float, float>> & filters_parsed, std::tuple<std::string, float, float> & cur_parsed){
  std::string field = std::get<0>(cur_parsed);
  float lower_bound = std::get<1>(cur_parsed);
  float upper_bound = std::get<2>(cur_parsed);
  
  for (size_t i=0; i<filters_parsed.size(); i++){
    if (std::get<0>(filters_parsed[i]) == field){
      std::get<1>(filters_parsed[i]) = std::max(std::get<1>(filters_parsed[i]), lower_bound);
      std::get<2>(filters_parsed[i]) = std::min(std::get<2>(filters_parsed[i]), upper_bound);
      return;
    }
  }
  filters_parsed.push_back(cur_parsed);
}

std::tuple<std::string, float, float> parse_filter(std::string & filter){
  std::string field;
  float lower_bound = -std::numeric_limits<float>::infinity();
  float upper_bound = std::numeric_limits<float>::infinity();
  
  std::vector<std::string> comparators = {"<=", ">=", "==", "<", ">"};
  for (size_t i=0; i<comparators.size(); i++){
    std::string comparator = comparators[i];
    size_t comparator_loc = filter.find(comparator);
    if(comparator_loc == std::string::npos){
      continue;
    }
    field = filter.substr(0, comparator_loc);
    std::string bound = filter.substr(comparator_loc+comparator.length(), filter.length());

    if (comparator == "<"){
      upper_bound = std::stof(bound)*(1 - std::numeric_limits<float>::epsilon());
    }else if (comparator == ">"){
      lower_bound = std::stof(bound)*(1 + std::numeric_limits<float>::epsilon());
    }else if (comparator == "<="){
      upper_bound = std::stof(bound);
    }else if (comparator == ">="){
      lower_bound = std::stof(bound);
    }else if (comparator == "=="){
      upper_bound = std::stof(bound);
      lower_bound = std::stof(bound);
    }
    break;
  }

  return std::make_tuple(field, lower_bound, upper_bound); 
}