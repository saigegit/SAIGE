#ifndef LDMAT_HPP
#define LDMAT_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>





void setGlobalVarsInCPP_LDmat(std::string t_impute_method,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,

                               double t_missing_cutoff,
                               double t_maxMAFLimit,

                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
                               unsigned int t_max_markers_region);




void LDmatRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           arma::mat & annoIndicatorMat,
                           std::string t_outputFile,
                           unsigned int t_n,           // sample size
                           bool t_isImputation,
                           std::vector<std::string> & annoStringVec,
                           std::string regionName);


void writeOutfile_single_LDmat(
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
                        std::vector<std::string> & infoVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & missingRateVec,
                        std::vector<uint32_t> & N_Vec,
                        std::string regionName);


bool openOutfile_single_LDmat(bool isappend);

void closeOutfile_single_LDmat();
#endif
