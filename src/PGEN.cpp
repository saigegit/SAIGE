#include <boost/algorithm/string.hpp>
#include "PGEN.hpp"

#include <string>
#include <iostream>
#include <fstream>

namespace PGEN {

static const double kGenoRDoublePairs[32] ALIGNV16 = PAIR_TABLE16(0.0, 1.0, 2.0, NA_REAL);

PgenClass::PgenClass(std::string pgenFile,
    std::string psamFile,
    std::string pvarFile,
    std::vector<std::string>& sampleInModel) : pgenFile(pgenFile), psamFile(psamFile), pvarFile(pvarFile), _allele_idx_offsetsp(nullptr)
{
    std::cout << "pgenFile" << pgenFile << std::endl;
    std::cout << "psamFile" << psamFile << std::endl;
    std::cout << "pvarFile" << pvarFile << std::endl;
    readPvarFile();
    readPsamFile();
    setPosSampleInPgen(sampleInModel);

    N0 = sampleInPgen.size();
    N = posSampleInPgen.size();
    M = variantId.size();

    loadPgen(N0);
    readBuffer.resize(_subset_size);
    
}

void PgenClass::loadPgen(uint32_t raw_sample_ct)
{
    // similar to pgenlibr/src/pgenlibr.cpp
    plink2::PgenHeaderCtrl header_ctrl;
    _info_ptr = std::unique_ptr<plink2::PgenFileInfo>(new plink2::PgenFileInfo);
    uintptr_t pgfi_alloc_cacheline_ct;
    char errstr_buf[plink2::kPglErrstrBufBlen];

    plink2::PreinitPgfi(_info_ptr.get());
    uint32_t cur_sample_ct = raw_sample_ct;
    uint32_t cur_variant_ct = UINT32_MAX;
    if(plink2::PgfiInitPhase1(pgenFile.c_str(), nullptr, cur_variant_ct, cur_sample_ct, &header_ctrl, _info_ptr.get(), &pgfi_alloc_cacheline_ct, errstr_buf) != plink2::kPglRetSuccess) {
        Rcpp::stop(&(errstr_buf[7]));
    }
    const uint32_t raw_variant_ct = _info_ptr->raw_variant_ct;
    if(header_ctrl & 0x30) {
        _allele_idx_offsetsp = plink2::CreateRefcountedWptr(cur_variant_ct + 1);
        _info_ptr->allele_idx_offsets = _allele_idx_offsetsp->p;
    }
    _info_ptr->max_allele_ct = 2;
    if((header_ctrl & 0xc0) == 0xc0) {
        const uintptr_t raw_variant_ctl = plink2::DivUp(cur_variant_ct, plink2::kBitsPerWord);
        _nonref_flagsp = plink2::CreateRefcountedWptr(raw_variant_ctl + 1);
        _info_ptr->nonref_flags = _nonref_flagsp->p;
    }
    const uint32_t file_sample_ct = _info_ptr->raw_sample_ct;
    unsigned char* pgfi_alloc = nullptr;
    if(plink2::cachealigned_malloc(pgfi_alloc_cacheline_ct * plink2::kCacheline, &pgfi_alloc)) {
        Rcpp::stop("Out of memory");
    }
    uint32_t max_vrec_width;
    uintptr_t pgr_alloc_cacheline_ct;
    if(plink2::PgfiInitPhase2(header_ctrl, 1, 0, 0, 0, raw_variant_ct, &max_vrec_width, _info_ptr.get(), pgfi_alloc, &pgr_alloc_cacheline_ct, errstr_buf)) {
        if(pgfi_alloc && (!_info_ptr->vrtypes)) {
            plink2::aligned_free(pgfi_alloc);
        }
        Rcpp::stop(&(errstr_buf[7]));
    }
    if((!_allele_idx_offsetsp) && (_info_ptr->gflags & 4)) {
        Rcpp::stop("Multiallelic variants and phase/dosage info simultaneously present; not supported by SAIGE");
    }
    _state_ptr = std::unique_ptr<plink2::PgenReader>(new plink2::PgenReader);

    if(! _state_ptr) {
        Rcpp::stop("Out of memory");
    }
    plink2::PreinitPgr(_state_ptr.get());
    plink2::PgrSetFreadBuf(nullptr, _state_ptr.get());
    const uintptr_t pgr_alloc_main_byte_ct = pgr_alloc_cacheline_ct * plink2::kCacheline;
    const uintptr_t sample_subset_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerVec) * plink2::kBytesPerVec;
    const uintptr_t cumulative_popcounts_byte_ct = plink2::DivUp(file_sample_ct, plink2::kBitsPerWord * plink2::kInt32PerVec) * plink2::kBytesPerVec;
    const uintptr_t genovec_byte_ct = plink2::DivUp(file_sample_ct, plink2::kNypsPerVec) * plink2::kBytesPerVec;
    const uint32_t is_not_plink1_bed = (_info_ptr->vrtypes != nullptr);
    uintptr_t raregeno_byte_ct = 0;
    uintptr_t difflist_sample_ids_byte_ct = 0;
    if(is_not_plink1_bed) {
        const uint32_t max_stored_single_difflist_len = file_sample_ct / plink2::kPglMaxDifflistLenDivisor;
        raregeno_byte_ct = plink2::DivUp(2 * max_stored_single_difflist_len, plink2::kNypsPerVec) * plink2::kBytesPerVec;
        difflist_sample_ids_byte_ct = plink2::RoundUpPow2(3 * max_stored_single_difflist_len * sizeof(int32_t), plink2::kBytesPerVec);
    }
    const uintptr_t ac_byte_ct = plink2::RoundUpPow2(file_sample_ct * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    const uintptr_t ac2_byte_ct = plink2::RoundUpPow2(file_sample_ct * 2 * sizeof(plink2::AlleleCode), plink2::kBytesPerVec);
    uintptr_t multiallelic_hc_byte_ct = 0;
    if(_info_ptr->max_allele_ct != 2) {
        multiallelic_hc_byte_ct = 2 * sample_subset_byte_ct + ac_byte_ct + ac2_byte_ct;
    }
    const uintptr_t dosage_main_byte_ct = plink2::DivUp(file_sample_ct, (2 * plink2::kInt32PerVec)) * plink2::kBytesPerVec;
    unsigned char* pgr_alloc;
    if(plink2::cachealigned_malloc(pgr_alloc_main_byte_ct + (2 * plink2::kPglNypTransposeBatch + 5) * sample_subset_byte_ct + cumulative_popcounts_byte_ct + (1 + plink2::kPglNypTransposeBatch) * genovec_byte_ct + raregeno_byte_ct + difflist_sample_ids_byte_ct + multiallelic_hc_byte_ct + dosage_main_byte_ct + plink2::kPglBitTransposeBufbytes + 4 * (plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8), &pgr_alloc)) {
        Rcpp::stop("Out of memory");
    }
    plink2::PglErr reterr = plink2::PgrInit(pvarFile.c_str(), max_vrec_width, _info_ptr.get(), _state_ptr.get(), pgr_alloc);
    if(reterr != plink2::kPglRetSuccess) {
        if(!plink2::PgrGetFreadBuf(_state_ptr.get())) {
            plink2::aligned_free(pgr_alloc);
        }
        snprintf(errstr_buf, plink2::kPglErrstrBufBlen, "PgrInit() error %d", static_cast<int>(reterr));
        Rcpp::stop(errstr_buf);
    }
    unsigned char* pgr_alloc_iter = &(pgr_alloc[pgr_alloc_main_byte_ct]);
    _subset_include_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _subset_include_interleaved_vec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);

#ifdef USE_AVX2
    _subset_include_interleaved_vec[-3] = 0;
    _subset_include_interleaved_vec[-2] = 0;
#endif
    _subset_include_interleaved_vec[-1] = 0;

    _subset_cumulative_popcounts = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[cumulative_popcounts_byte_ct]);
    _pgv.genovec = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[genovec_byte_ct]);
    if(is_not_plink1_bed) {
        _raregeno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[raregeno_byte_ct]);
        _difflist_sample_ids_buf = reinterpret_cast<uint32_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[difflist_sample_ids_byte_ct]);
    } else {
        _raregeno_buf = nullptr;
        _difflist_sample_ids_buf = nullptr;
    }
    if(multiallelic_hc_byte_ct) {
        _pgv.patch_01_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_01_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac_byte_ct]);
        _pgv.patch_10_set = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
        _pgv.patch_10_vals = reinterpret_cast<plink2::AlleleCode*>(pgr_alloc_iter);
        pgr_alloc_iter = &(pgr_alloc_iter[ac2_byte_ct]);
    } else {
        _pgv.patch_01_set = nullptr;
        _pgv.patch_01_vals = nullptr;
        _pgv.patch_10_set = nullptr;
        _pgv.patch_10_vals = nullptr;
    }
    _pgv.phasepresent = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.phaseinfo = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_present = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[sample_subset_byte_ct]);
    _pgv.dosage_main = reinterpret_cast<uint16_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[dosage_main_byte_ct]);
    _subset_size = file_sample_ct;

    _transpose_batch_buf = reinterpret_cast<plink2::VecW*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglBitTransposeBufbytes]);
    _multivar_vmaj_geno_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * genovec_byte_ct]);
    _multivar_vmaj_phasepresent_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_vmaj_phaseinfo_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * sample_subset_byte_ct]);
    _multivar_smaj_geno_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 4]);
    _multivar_smaj_phaseinfo_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
    _multivar_smaj_phasepresent_batch_buf = reinterpret_cast<uintptr_t*>(pgr_alloc_iter);
    // pgr_alloc_iter = &(pgr_alloc_iter[plink2::kPglNypTransposeBatch * plink2::kPglNypTransposeBatch / 8]);
}

void PgenClass::Read(std::vector<double>& buf, int variant_idx, int allele_idx) {
    // similar to pgenlibr/src/pgenlibr.cpp
    if(!_info_ptr.get()) {
        Rcpp::stop("pgen is closed");
    }
    if(static_cast<uint32_t>(variant_idx) >= _info_ptr->raw_variant_ct) {
        Rcpp::stop("variant_num out of range ("+std::to_string(variant_idx + 1)+"; must be 1.."+std::to_string(_info_ptr->raw_variant_ct)+")");
    }
    if(buf.size() != _subset_size) {
        Rcpp::stop("buf has wrong length ("+std::to_string(buf.size())+"; "+std::to_string(_subset_size)+" expected)");
    }
    uint32_t dosage_ct;
    plink2::PglErr reterr = plink2::PgrGet1D(_subset_include_vec, _subset_index, _subset_size, variant_idx, allele_idx, _state_ptr.get(), _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, &dosage_ct);
    if(reterr != plink2::kPglRetSuccess) {
        Rcpp::stop("PgrGet1D() error " + std::to_string(static_cast<int>(reterr)));
    }
    plink2::Dosage16ToDoubles(kGenoRDoublePairs, _pgv.genovec, _pgv.dosage_present, _pgv.dosage_main, _subset_size, dosage_ct, &buf[0]);
}


void to_upper(std::string str) {
    std::transform(str.begin(), str.end(), str.begin(), toupper);
}

void PgenClass::readPvarFile()
{
    std::cout << "Reading pvar file...." << std::endl;
    std::ifstream bim(pvarFile);
    std::string line;

    int col_pos = -1;
    int col_id = -1;
    int col_ref = -1;
    int col_alt = -1;

    while(getline(bim, line)) {
        std::vector<std::string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        boost::replace_all(line_elements.back(), "\r", "");

        if(line[0] == '#') {
            if(line.rfind("#CHROM", 0) == 0) {
                for(int i=0; i<line_elements.size(); i++) {
                    if(line_elements[i] == "POS") col_pos = i;
                    if(line_elements[i] == "ID") col_id = i;
                    if(line_elements[i] == "REF") col_ref = i;
                    if(line_elements[i] == "ALT") col_alt = i;
                }
                if(col_pos == -1 || col_id == -1 || col_ref == -1 || col_alt == -1) {
                    Rcpp::stop("pvar file does not have pos/id/ref/alt columns");
                }
            }
            continue;
        }

        chr.push_back(line_elements[0]);
        variantId.push_back(line_elements[col_id]);
        position.push_back(std::stoi(line_elements[col_pos]));
        to_upper(line_elements[col_ref]);
        to_upper(line_elements[col_alt]);
        ref.push_back(line_elements[col_ref]);
        alt.push_back(line_elements[col_alt]);
    }
}

void PgenClass::readPsamFile()
{
    std::cout << "Reading psam file (" << psamFile << ")...." << std::endl;
    std::ifstream fam(psamFile);

    std::string line;
    size_t idx = 0;
    while(getline(fam, line)) {
        if(line[0] == '#') continue;
        std::vector<std::string> line_elements;
        boost::split(line_elements, line, boost::is_any_of("\t "));
        sampleInPgen[line_elements[0]] = idx;  // IID
        idx++;
    }
}


void PgenClass::setPosSampleInPgen(std::vector<std::string>& sampleInModel)
{
    std::cout << "Setting position of samples in PGEN file...." << std::endl;
    auto m_N = sampleInModel.size();
    posSampleInPgen.reserve(m_N);

    for(auto& sample : sampleInModel) {
        auto res = sampleInPgen.find(sample);
        if(res == sampleInPgen.end()) {
            Rcpp::stop("At least one subject requested is not in pvar file("+sample+").");
        }
        posSampleInPgen.push_back(res->second);
    }    
}

void PgenClass::closegenofile() {
    if(_info_ptr) {
        CondReleaseRefcountedWptr(&_allele_idx_offsetsp);
        CondReleaseRefcountedWptr(&_nonref_flagsp);
        if(_info_ptr->vrtypes) {
            plink2::aligned_free(_info_ptr->vrtypes);
        }
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgfi(_info_ptr.get(), &reterr);
        _info_ptr.release();
    }
    if(_state_ptr) {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgr(_state_ptr.get(), &reterr);
        if(PgrGetFreadBuf(_state_ptr.get())) {
            plink2::aligned_free(PgrGetFreadBuf(_state_ptr.get()));
        }
        _state_ptr.release();
    }
    _subset_size = 0;
}


void PgenClass::getOneMarker(
    uint64_t& t_gIndex,        // different meanings for different genoType
    std::string& t_ref,       // REF allele
    std::string& t_alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
    std::string& t_marker,    // marker ID extracted from genotype file
    uint32_t& t_pd,           // base position
    std::string& t_chr,       // chromosome
    double& t_altFreq,        // frequency of ALT allele
    double& t_altCounts,      // counts of ALT allele
    double& t_missingRate,    // missing rate
    double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PGEN)
    bool& t_isOutputIndexForMissing,               // if true, output index of missing genotype data
    std::vector<uint>& t_indexForMissing,     // index of missing genotype data
    bool& t_isOnlyOutputNonZero,                   // if true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
    std::vector<uint>& t_indexForNonZero,     // only used when t_isOnlyOutputNonZero = TRUE
    arma::vec& OneMarkerG1)
{
    uint32_t numMissing = 0;

    t_indexForMissing.clear();
    t_indexForNonZero.clear();

    t_chr = chr[t_gIndex];
    t_pd = position[t_gIndex];
    t_ref = ref[t_gIndex];
    t_alt = alt[t_gIndex];
    t_marker = variantId[t_gIndex];
    
    // fill read_buffer with the number of ALT reads per individual
    int alleleIdx = 1;
    Read(readBuffer, t_gIndex, alleleIdx);

    uint j = 0;
    t_altCounts = 0;
    for(uint32_t i = 0; i < N; i++)
    {
        auto sampleIdx = posSampleInPgen[i];
        auto numAltReads = readBuffer[sampleIdx];
        if(std::isnan(numAltReads)) {
            numMissing++;
            if(t_isOutputIndexForMissing) {
                t_indexForMissing.push_back(i);
            }
        } else {
            t_altCounts += numAltReads;
        }

        if(numAltReads > 0) {
            t_indexForNonZero.push_back(i);
        }
        if(t_isOnlyOutputNonZero){
            if(numAltReads > 0) {
                OneMarkerG1[j] = numAltReads;
                j = j + 1;
            }
        }else{
            OneMarkerG1[i] = numAltReads;
        }
    }
    int count = N - numMissing;
    t_missingRate = (double)numMissing / (double)N;
    t_imputeInfo = 1;

    if(count > 0){
            t_altFreq = t_altCounts / (double)count / 2;
    }else{
            t_altFreq = 0;
    }

    if(t_isOnlyOutputNonZero){
        OneMarkerG1.resize(j);
    }
}

}