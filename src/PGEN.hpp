#pragma once

// All C++ codes related to PLINK file manipulation

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <unordered_map>
#include <memory>
#include <pgenlib_read.h>
#include <pvar_ffi_support.h>
#include <pgenlib_ffi_support.h>

namespace PGEN {

typedef unsigned int uint;

class PgenClass {
private:
    // information from pvar file
    std::vector<std::string> chr;              // Chromosome code (either an integer, or 'X'/'Y'/'XY'/'MT'; '0' indicates unknown) or name
    std::vector<std::string> variantId;        // Variant identifier
    std::vector<uint32_t> position;            // Base-pair coordinate (1-based; limited to 2^31-2)
    std::vector<std::string> ref;              // Allele 1
    std::vector<std::string> alt;              // Allele 2
    size_t M;

    // information from psam file
    std::unordered_map<std::string, size_t> sampleInPgen;
    uint32_t N, N0;

    // PGEN file
    std::unique_ptr<plink2::PgenFileInfo> _info_ptr;
    std::unique_ptr<plink2::PgenReader> _state_ptr;
    std::string pgenFile, psamFile, pvarFile;
    std::vector<size_t> posSampleInPgen;

    plink2::RefcountedWptr* _allele_idx_offsetsp;
    plink2::RefcountedWptr* _nonref_flagsp;
    uintptr_t* _subset_include_vec;
    uintptr_t* _subset_include_interleaved_vec;
    uint32_t* _subset_cumulative_popcounts;
    plink2::PgrSampleSubsetIndex _subset_index;
    uint32_t _subset_size;
    plink2::PgenVariant _pgv;
    uintptr_t* _raregeno_buf;
    uint32_t* _difflist_sample_ids_buf;
    plink2::VecW* _transpose_batch_buf;
    uintptr_t* _multivar_vmaj_geno_buf;
    uintptr_t* _multivar_vmaj_phasepresent_buf;
    uintptr_t* _multivar_vmaj_phaseinfo_buf;
    uintptr_t* _multivar_smaj_geno_batch_buf;
    uintptr_t* _multivar_smaj_phaseinfo_batch_buf;
    uintptr_t* _multivar_smaj_phasepresent_batch_buf;
    std::vector<double> readBuffer;

    // internal functions
    void readPvarFile();
    void readPsamFile();
    void loadPgen(uint32_t raw_sample_ct);
public:

    PgenClass(std::string pgenFile,
        std::string psamFile,
        std::string pvarFile,
        std::vector<std::string>& sampleInModel);


    void setPosSampleInPgen(std::vector<std::string>& SampleInModel);
    std::vector<uint32_t> getPosMarkerInPgen(std::vector<std::string> MarkerReqstd);
    void Read(std::vector<double>& buf, int variant_idx, int allele_idx);

    void getOneMarker(
        uint64_t& t_gIndex,        // different meanings for different genoType
        std::string& t_ref,       // REF allele
        std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
        std::string& t_marker,    // marker ID extracted from genotype file
        uint32_t& t_pd,           // base position
        std::string& t_chr,       // chromosome
        double& t_altFreq,        // frequency of ALT allele
        double& t_altCounts,      // counts of ALT allele
        double& t_missingRate,    // missing rate
        double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
        bool& t_isOutputIndexForMissing,               // if true, output index of missing genotype data
        std::vector<uint>& t_indexForMissing,     // index of missing genotype data
        bool& t_isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
        std::vector<uint>& t_indexForNonZero,
        arma::vec& OneMarkerG1);
    void getOneMarker(uint64_t t_gIndex,
        double& t_altFreq,
        double& t_missingRate,
        std::string& t_chr,
        arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForMissing, indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        getOneMarker(t_gIndex, ref, alt, marker, pd, t_chr, t_altFreq, altCounts, t_missingRate, imputeInfo,
            isOutputIndexForMissing, indexForMissing, isOnlyOutputNonZero, indexForNonZero, OneMarkerG1);
    }

    void getOneMarker(uint64_t t_gIndex,
        double& t_altFreq,
        double& t_missingRate,
        std::vector<uint>& t_indexForMissing,
        arma::vec& OneMarkerG1)
    {
        std::string ref, alt, marker, chr;
        uint32_t pd;
        double altCounts, imputeInfo;
        std::vector<uint> indexForNonZero;
        bool isOutputIndexForMissing = false;
        bool isOnlyOutputNonZero = false;
        getOneMarker(t_gIndex, ref, alt, marker, pd, chr, t_altFreq, altCounts, t_missingRate, imputeInfo,
            isOutputIndexForMissing, t_indexForMissing, isOnlyOutputNonZero, indexForNonZero, OneMarkerG1);
    }


    uint32_t getN0() { return N0; }
    uint32_t getN() { return N; }
    uint32_t getM() { return M; }


    void closegenofile();

};

}