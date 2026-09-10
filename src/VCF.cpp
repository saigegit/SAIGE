// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
 
//#include "variant_group_iterator.hpp"

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

using namespace std;

 
namespace VCF {

void stride_reduce_dosage(int32_t *values, int n_samples, size_t max_ploidy) {
    // similar to the example in htslib/vcf.h
    for(size_t i=0; i<n_samples; i++) {
        int32_t *ptr = values + i*max_ploidy;
        bool missing = true;
        int num_alt_alleles = 0;
        for (int j=0; j<max_ploidy; j++) {
            // if true, the sample has smaller ploidy
            if(ptr[j]==bcf_int32_vector_end) break;
            // missing allele
            if(bcf_gt_is_missing(ptr[j])) continue;
            missing = false;
            num_alt_alleles += (bcf_gt_allele(ptr[j]) > 0);
        }
        if(missing) {
            values[i] = -1;
        } else {
            values[i] = num_alt_alleles;
        }
    }
}

VcfClass::VcfClass(std::string vcfFileName,
                   std::string fmtField,
                   std::string vcfFilterString,
                   bool isSparseDosageInVcf,
                   std::vector<std::string> SampleInModel) : vcfFileName(vcfFileName), fmtField(fmtField), vcf_record(nullptr), vcf_iterator(nullptr), has_line(0), iterate_over_all_lines(true), use_sites(false), vcf_idx(0), vcf_filter_string(vcfFilterString), multiallelic_vcf_record(nullptr), multiallelic_idx(-1) {
    vcf_string = {0,0,0};
    vcf_record = bcf_init();

    // filename could contain a range: "exome/vcf/[range:0000000012-0000000025].vcf.bgz"
    range_pos = vcfFileName.find("[range:");
    multiple_vcf_files = (range_pos != string::npos);
    if(multiple_vcf_files) {
        multiple_vcf_files = true;
        auto sep_pos = vcfFileName.find("-", range_pos);
        auto end_pos = vcfFileName.find("]", range_pos);
        auto range_start = vcfFileName.substr(range_pos + 7, sep_pos - range_pos - 7);  // 0000000012
        auto range_end = vcfFileName.substr(sep_pos + 1, end_pos - sep_pos - 1);        // 0000000025
        range_len = range_start.size();
        vcfFileName_first = vcfFileName.substr(0, range_pos) + range_start + vcfFileName.substr(end_pos + 1);
        vcfFileName_last  = vcfFileName.substr(0, range_pos) + range_end   + vcfFileName.substr(end_pos + 1);
        VcfClass::vcfFileName = vcfFileName_first;
    }

    vcf_file.resize(1);
    vcf_index.resize(1);
    vcf_header.resize(1);
    vcf_filter.resize(1);
    auto isVcfOpen = setVcfObj(0);
    if(!isVcfOpen) {
        std::cout << "Failed opening VCF file" << std::endl;
        return;
    }
    setPosSampleInVcf(SampleInModel);

   
    vcfnorm_args.tmp_lines = nullptr;
    vcfnorm_args.lines = nullptr;
    vcfnorm_args.ntmp_lines = 0;
    vcfnorm_args.mtmp_lines = 0;
    vcfnorm_args.als = nullptr;
    vcfnorm_args.tmp_arr1 = nullptr;
    vcfnorm_args.ntmp_arr1 = 0;
    vcfnorm_args.hdr = nullptr;
    vcfnorm_args.out_hdr = nullptr;
    vcfnorm_args.force = 1;
    vcfnorm_args.force_warned = 1;
    vcfnorm_args.keep_sum_ad = 0;
    vcfnorm_args.ma_use_ref_allele = 0;
    
    bcftools::init_data(&vcfnorm_args);
}
VcfClass::~VcfClass() {
    if(vcf_string.m) {
        free(vcf_string.s);
    }
    if(multiallelic_vcf_record != nullptr) {
        bcf_destroy(multiallelic_vcf_record);
    } else {
        bcf_destroy(vcf_record);
    }
    for(size_t i=0; i<vcf_file.size(); i++) {
        closeAll(i);
    }
    if(vcfnorm_args.hdr != nullptr) {
        bcf_hdr_destroy(vcfnorm_args.hdr);
    }
    bcftools::destroy_data(&vcfnorm_args);
}

void VcfClass::closeAll(size_t idx) {
    if(vcf_header[idx] != nullptr) {
        bcf_hdr_destroy(vcf_header[idx]);
        vcf_header[idx] = nullptr;
    }
    if(vcf_file[idx] != nullptr) {
        hts_close(vcf_file[idx]);
        vcf_file[idx] = nullptr;
    }
    if(vcf_index[idx] != nullptr) {
        tbx_destroy(vcf_index[idx]);
        vcf_index[idx] = nullptr;
    }
    if(vcf_iterator != nullptr) {
        bcf_itr_destroy(vcf_iterator);
        vcf_iterator = nullptr;
    }
    if(vcf_filter[idx] != nullptr) {
        filter_destroy(vcf_filter[idx]);
        vcf_filter[idx] = nullptr;
    }
}

bool VcfClass::setVcfObj(size_t idx) {
    closeAll(idx);

    vcf_file[idx] = vcf_open(vcfFileName.c_str(), "r");
    if(vcf_file[idx] == nullptr) {
        std::cerr << "WARNING: Open VCF " << vcfFileName << " failed" << std::endl;
        return false;
    }
    
    vcf_header[idx] = bcf_hdr_read(vcf_file[idx]);
    if(vcf_header[idx] == nullptr) {
        std::cerr << "WARNING: Reading VCF header of " << vcfFileName << " failed" << std::endl;
        hts_close(vcf_file[idx]);
        return false;
    }

    vcf_index[idx] = tbx_index_load(vcfFileName.c_str());

    std::cout << "Open VCF " << vcfFileName << " done" << std::endl;
    cout << "To read the field " << fmtField << endl;
    std::cout << "Number of meta lines in the vcf file (lines starting with ##): " << vcf_header[idx]->nhrec << endl;
    N0 = bcf_hdr_nsamples(vcf_header[idx]);
    std::cout << "Number of samples in the vcf file: " << N0 << endl;

    if(bcf_hdr_id2int(vcf_header[idx], BCF_DT_ID, fmtField.c_str()) >= 0) {
        // Regarding dosage fields:
        // DS and GT are equivalent.
        // That is, a given VCF may list its dosage field as "GT" rather than "DS".
        if((fmtField == "DS" || fmtField == "GT") && bcf_hdr_id2int(vcf_header[idx], BCF_DT_ID, "HDS") >= 0) {
            fmtField = "HDS";
        }
    } else {
        std::cerr << "ERROR: fmtField (" << fmtField << ") not present in genotype file." << std::endl;
        return false;
    }

    if(vcf_filter_string != ""){
        vcf_filter[idx] = bcftools::filter_init(vcf_header[idx], vcf_filter_string.c_str());
    }
    
    return true;
}

void VcfClass::set_iterator(std::string& variantList) {
    // determine number of sites
    size_t num_sites = 0;
    for(size_t i=0; i<variantList.size(); i++) num_sites += (variantList[i] == '\t');

    sites.clear();
    istringstream iss(variantList);
    string site_name;
    getline(iss, site_name, '\t');
    for(size_t i=0; i<num_sites; i++) {
        site_t site;
        string str_pos;
        getline(iss, site.chrom, ':');
        getline(iss, str_pos, ':');
        site.pos = std::stoi(str_pos);
        getline(iss, site.ref, ':');
        getline(iss, site.alt, '\t');
        sites.insert(site);
    }

    // open all vcf files to efficiently iterate through sites at a (modest) memory increase
    if(multiple_vcf_files && vcf_file.size() == 1) {
        size_t num_files = 1 + stoi(vcfFileName_last.substr(range_pos, range_len)) - stoi(vcfFileName_first.substr(range_pos, range_len));
        vcf_file.resize(num_files);
        vcf_index.resize(num_files);
        vcf_header.resize(num_files);
        vcf_filter.resize(num_files);
        vcf_first_line.resize(num_files);
        for(size_t i=1; i<num_files; i++) {
            nextVcfFileName();
            auto isVcfOpen = setVcfObj(i);
            if(!isVcfOpen) {
                std::cout << "Failed opening VCF file" << std::endl;
                return;
            }
            setPosSampleInVcf(i);
        }
        for(size_t i=0; i<num_files; i++) {
            auto ret_val = bcf_read(vcf_file[i], vcf_header[i], vcf_record);
            if(ret_val == 0) {
                bcf_unpack(vcf_record, BCF_UN_STR);
                vcf_first_line[i] = { bcf_hdr_id2name(vcf_header[i], vcf_record->rid),
                                      vcf_record->pos + 1,
                                      vcf_record->d.allele[0],
                                      vcf_record->d.allele[1] };
            }
        }
        vcfFileName = vcfFileName_first;
    }

    has_line = true;
    iterate_over_all_lines = false;
    use_sites = true;
    vcf_idx = 0;
    set_iterator();
    move_forward_iterator(1);
}

void VcfClass::set_iterator(std::string& chrom, int& beg_pd, int& end_pd) {
    if(chrom == "entire_file") {
        iterator_string = "not used";
        iterate_over_all_lines = true;
    } else {
        if(vcf_index[vcf_idx] == nullptr) {
            Rcpp::stop("ERROR: Open VCF index for " + vcfFileName + " failed and not processing entire file");
        }
        iterator_string = chrom + ":" + to_string(beg_pd) + "-" + to_string(end_pd);
        iterate_over_all_lines = false;
    }
    has_line = true;
    set_iterator();
    move_forward_iterator(1);
}

void VcfClass::set_iterator() {
    if(! iterate_over_all_lines) {
        if(vcf_iterator != nullptr) {
            bcf_itr_destroy(vcf_iterator);
        }
        if(use_sites) {
            vector<string> sites_locus;
            vector<char*> sites_c_array;
            sites_locus.reserve(sites.size());
            sites_c_array.reserve(sites.size());
            for (const auto& site : sites) {
                sites_locus.push_back(site.chrom + ":" + to_string(site.pos) + "-" + to_string(site.pos));
                sites_c_array.push_back(const_cast<char*>(sites_locus.back().c_str()));
            }
            vcf_iterator = tbx_itr_regarray(vcf_index[vcf_idx], &sites_c_array[0], sites.size());
        } else {
            vcf_iterator = tbx_itr_querys(vcf_index[vcf_idx], iterator_string.c_str());
        }

        if(vcf_iterator == nullptr) {
            std::cerr << "ERROR: invalid region " << iterator_string << std::endl;
        }
    }
}

void VcfClass::move_forward_iterator(int i) {
    if(! has_line) return;
    for(int j = 0; j < i; j++){
        while(true) {
            if(multiallelic_vcf_record != nullptr && multiallelic_idx == multiallelic_vcf_record->n_allele - 2) {
                vcf_record = multiallelic_vcf_record;
                vcf_header[vcf_idx] = vcfnorm_args.hdr;
                multiallelic_vcf_record = nullptr;
                bcf_hdr_destroy(vcfnorm_args.out_hdr);
                vcfnorm_args.out_hdr = nullptr;
            }
            if(multiallelic_vcf_record != nullptr) {
                multiallelic_idx++;
                vcf_record = vcfnorm_args.tmp_lines[multiallelic_idx];
            } else {
                move_forward_iterator();
                if(! has_line) return;
                bcf_unpack(vcf_record, BCF_UN_STR);
                if(vcf_record->n_allele > 2) {
                    bcf_unpack(vcf_record, BCF_UN_ALL);
                    multiallelic_vcf_record = vcf_record;
                    vcfnorm_args.hdr = vcf_header[vcf_idx];
                    vcfnorm_args.out_hdr = bcf_hdr_dup(vcf_header[vcf_idx]);
                    vcf_header[vcf_idx] = vcfnorm_args.out_hdr;
                    bcftools::split_multiallelic_to_biallelics(&vcfnorm_args, multiallelic_vcf_record);
                    multiallelic_idx = 0;
                    vcf_record = vcfnorm_args.tmp_lines[multiallelic_idx];
                }
            }
            if(! use_sites) {
                break;
            }
            string chrom = bcf_hdr_id2name(vcf_header[vcf_idx], vcf_record->rid);
            if(sites.count({chrom, vcf_record->pos + 1, vcf_record->d.allele[0], vcf_record->d.allele[1]})) {
                break;
            }
        }
    }
    bcf_unpack(vcf_record, BCF_UN_ALL);
}


void VcfClass::move_forward_iterator() {
    int ret_val;
    if(iterate_over_all_lines) {
        ret_val = bcf_read(vcf_file[vcf_idx], vcf_header[vcf_idx], vcf_record);
    } else {
        ret_val = tbx_itr_next(vcf_file[vcf_idx], vcf_index[vcf_idx], vcf_iterator, &vcf_string); 
        vcf_parse(&vcf_string, vcf_header[vcf_idx], vcf_record);
    }
    
    // ret_val is 0 on success; -1 on end of file; < -1 on critical error (for both bcf_read and tbx_itr_next)
    has_line = (ret_val >= 0);
    while(! has_line) {
        if(ret_val < -1) {
            std::cout << "htslib error: " << vcf_record->errcode << endl;
            Rcpp::stop("Error reading VCF file.");
        }
        if(!multiple_vcf_files || vcfFileName == vcfFileName_last) {
            std::cout << "Reached the end of the VCF file" << std::endl;
            return;
        }
        if(use_sites) {
            vcf_idx++;
            bool file_contains_site = false;
            while(vcf_idx < vcf_file.size()) {
                for(auto& site : sites) {
                    file_contains_site = site.chrom == vcf_first_line[vcf_idx].chrom
                                    && site.pos >= vcf_first_line[vcf_idx].pos
                                    && (vcf_idx == vcf_file.size() - 1 || site.pos <= vcf_first_line[vcf_idx+1].pos);
                    if(file_contains_site) break;
                }
                if(file_contains_site) break;
                vcf_idx++;
            }
            if(! file_contains_site) {
                has_line = false;
                return;
            }
            set_iterator();
        } else {
            nextVcfFileName();
            auto isVcfOpen = setVcfObj(0);
            if(!isVcfOpen) {
                std::cout << "Failed opening VCF file" << std::endl;
                return;
            }
            setPosSampleInVcf(0);
            set_iterator();
        }
        
        if(iterate_over_all_lines) {
            ret_val = bcf_read(vcf_file[vcf_idx], vcf_header[vcf_idx], vcf_record);
        } else {
            ret_val = tbx_itr_next(vcf_file[vcf_idx], vcf_index[vcf_idx], vcf_iterator, &vcf_string);
            vcf_parse(&vcf_string, vcf_header[vcf_idx], vcf_record);
        }
        has_line = (ret_val >= 0);
    }

}

bool VcfClass::check_iterator_end() {
    return !has_line;
}

void VcfClass::setPosSampleInVcf(size_t vcf_idx) {
    int num_samples = SampleInModel.size();
    string samples;
    for(size_t i=0; i<num_samples; i++) {
        if(samples.size() > 0) samples += ",";
        samples += SampleInModel[i];
    }

    int ret = bcf_hdr_set_samples(vcf_header[vcf_idx], samples.c_str(), 0);
    if(vcf_filter_string != ""){
        filter_destroy(vcf_filter[vcf_idx]);
        vcf_filter[vcf_idx] = bcftools::filter_init(vcf_header[vcf_idx], vcf_filter_string.c_str());
    }

    if(ret < 0) {
        Rcpp::stop("Error setting VCF samples");
    }
    if(ret > 0) {
        Rcpp::stop("At least one subject requested is not in VCF file.");
    }

    if(vcf_idx > 0) {
        std::vector<std::string> SampleInVcf = SampleInVcf;
        getSampleIDlist_vcfMatrix(vcf_idx);
        if(SampleInVcf != SampleInVcf) {
            Rcpp::stop("Processing multiple VCF files with a sample list requires that all VCF files have the same samples in the same order.");
        }
        return;
    }

    getSampleIDlist_vcfMatrix(0);
    if(N == 0) {
        SampleInModel = SampleInVcf;
    }

    posSampleInModel.reserve(N0);
    map<string, size_t> sample_to_model_idx;
    for(size_t i=0; i<N; i++) {
        sample_to_model_idx[SampleInModel[i]] = i;
    }
    
    for(auto& sample : SampleInVcf) {
        auto it = sample_to_model_idx.find(sample);
        if(it == sample_to_model_idx.end()) {
            posSampleInModel.push_back(-1);
        } else {
            posSampleInModel.push_back(it->second);
        }
    }
}

void VcfClass::setPosSampleInVcf(std::vector<std::string>& SampleInModel) {
    std::cout << "Setting position of samples in VCF files...." << std::endl;

    N = SampleInModel.size();
    VcfClass::SampleInModel = SampleInModel;
    setPosSampleInVcf(0);
}

// get dosages/genotypes of one marker   
void VcfClass::getOneMarker(
                        std::string& ref,       // REF allele
                        std::string& alt,       // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                        std::string& marker,    // marker ID extracted from genotype file
                        uint32_t& pd,           // base position
                        std::string& chr,       // chromosome
                        double& altFreq,        // frequency of ALT allele
                        double& altCounts,      // counts of ALT allele
                        double& missingRate,    // missing rate
                        double& imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                        bool isOutputIndexForMissing,               // if true, output index of missing genotype data
                        std::vector<uint32_t>& indexForMissingforOneMarker,     // index of missing genotype data
                        bool isOnlyOutputNonZero,                   // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                        std::vector<uint32_t>& indexForNonZero,
                        //if true, the marker has been read successfully
                        bool& isBoolRead,
                        arma::vec& dosages,
                        bool isImputation) {
    isBoolRead = true;
    if(! has_line) {
        isBoolRead = false;
        return;
    }

    if (vcf_record->errcode) {
        isBoolRead = false; 
        cout << "Failed reading VCF line, error " << vcf_record->errcode << std::endl;
        return;
    }

    chr = bcf_hdr_id2name(vcf_header[vcf_idx], vcf_record->rid);
    pd = vcf_record->pos + 1;
    ref = vcf_record->d.allele[0];
    alt = vcf_record->d.allele[1];
    for(size_t i=2; i<vcf_record->n_allele; i++) {
        alt += std::string(",") + vcf_record->d.allele[i];
    }
    //if(pd > 180020) exit(0);
    marker = vcf_record->d.id;

    if(isImputation) {
        float* r2;
        int r2_size;
        int ret = bcf_get_info_float(vcf_header[vcf_idx], vcf_record, "R2", &r2, &r2_size);
        if(r2_size == 0) {
            Rcpp::stop("Cannot impute since R2 field is missing.");
        }
        imputeInfo = r2[0];
        free(r2);
    } else {
        imputeInfo = 1.0;
    }
    
    double dosage;
    altCounts = 0;
    int missing_cnt = 0;
    
    dosages.clear();
    dosages.set_size(N);
    dosages.fill(arma::fill::zeros);
    
    int32_t *variant_dosages = NULL;
    int variant_dosages_size = 0;
    auto n_variant_dosages = bcf_get_genotypes(vcf_header[vcf_idx], vcf_record, &variant_dosages, &variant_dosages_size);
    std::size_t ploidy = N ? n_variant_dosages / N : 1;
    stride_reduce_dosage(variant_dosages, N, ploidy);

    const uint8_t* filter_samples;
    if(vcf_filter_string != "") {
        bcftools::filter_test(vcf_filter[vcf_idx], vcf_record, &filter_samples);
    }
    
    for(size_t i=0; i<N; i++) {
        if(vcf_filter_string != "" && filter_samples[i]) {
            variant_dosages[i] = -1;
        }

        dosages[posSampleInModel[i]] = variant_dosages[i];
        if(variant_dosages[i] == -1) {
            ++missing_cnt;
            indexForMissingforOneMarker.push_back(posSampleInModel[i]);
        } else if(variant_dosages[i] > 0) {
            indexForNonZero.push_back(posSampleInModel[i]);
            altCounts += variant_dosages[i];
        }
    }

    if(missing_cnt > 0) {
        if(missing_cnt == N) {
            altFreq = 0;
        } else {
            altFreq = altCounts / 2 / (double)(N - missing_cnt);
        }
        missingRate = missing_cnt/double(N);
    } else {
        altFreq = altCounts / 2 / (double)(N);
        missingRate = 0; 
    }
    free(variant_dosages);
}
 
void VcfClass::setIsSparseDosageInVcf(bool isSparseDosageInVcf){
    // not used
}

void VcfClass::getSampleIDlist_vcfMatrix(size_t idx){
    SampleInVcf = std::vector<std::string>(vcf_header[idx]->samples, vcf_header[idx]->samples + bcf_hdr_nsamples(vcf_header[idx]));
}

void VcfClass::nextVcfFileName() {
    size_t pos = range_pos + range_len - 1;
    vcfFileName[pos]++;
    while(vcfFileName[pos] > '9' && pos > 0) {
        vcfFileName[pos] = '0';
        pos--;
        vcfFileName[pos]++;
    }
}

}