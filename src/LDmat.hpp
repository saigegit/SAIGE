#ifndef LDMAT_HPP
#define LDMAT_HPP

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// HighFive: C++ wrapper for HDF5
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5PropertyList.hpp> 

using namespace HighFive;

// Declare your global HDF5 file handle if shared
// HighFive::File h5File;  // defined in LDmat.cpp

// === Existing declarations ===

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
        std::string t_genoType,
        std::vector<std::string> & t_genoIndex_prev,
        std::vector<std::string> & t_genoIndex,
        arma::mat & annoIndicatorMat,
        std::string t_outputFile,
        unsigned int t_n,
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


// === New HDF5-related functions ===

SEXP initHDF5(std::string h5Filename);

void writeSparseToHDF5(SEXP file_xptr,
                       std::string regionName,
                       std::vector<unsigned int>& rowVec,
                       std::vector<unsigned int>& colVec,
                       std::vector<int>& valVec);
#endif  // LDMAT_HPP
