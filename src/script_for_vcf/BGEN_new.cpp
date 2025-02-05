// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]

#include <cstring>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <set>
#include <numeric>
#include <utility>
#include <stdexcept>
#include <memory>
#include <time.h>
#include <stdint.h>
#include <zlib.h>

#include <cassert>
#include "../thirdParty/bgen/genfile/include/genfile/bgen/bgen.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/View.hpp"
#include "../thirdParty/bgen/genfile/include/genfile/bgen/IndexQuery.hpp"

#include "UTIL.hpp"
// #include <boost/iostreams/filter/zstd.hpp>
// #include "zstd.h"
// #include <boost/date_time.hpp>

#include <Rcpp.h>
#include "BGEN.hpp"

namespace BGEN {

BgenClass::BgenClass(std::string t_bgenFileName,
                     std::string t_bgenFileIndex,
                     std::vector<std::string>  t_SampleInBgen,
                     std::vector<std::string>  t_SampleInModel,
                     bool t_isSparseDosageInBgen,
                     bool t_isDropmissingdosagesInBgen,
                     std::string t_AlleleOrder, 
		     std::vector< std::string > const& t_ids_to_include,
		     Rcpp::DataFrame & t_ranges_to_include)      // added by Wenjian Bi on 03/14/2021: "ref-first" or "alt-first"
{
  setBgenObj(t_bgenFileName,
             t_bgenFileIndex,
             t_SampleInBgen,
	     t_ids_to_include,
	     t_ranges_to_include);
  
  setPosSampleInBgen(t_SampleInModel);
  
  setIsDropMissingDosagesInBgen(t_isDropmissingdosagesInBgen);
  
  setIsSparseDosageInBgen(t_isSparseDosageInBgen);
  m_AlleleOrder = t_AlleleOrder;
}


void BgenClass::setBgenObj(const std::string t_bgenFileName,
                           const std::string t_bgenFileIndex,
                           std::vector<std::string> & t_SampleInBgen,
			   std::vector< std::string > const& t_ids_to_include,
			   Rcpp::DataFrame & t_ranges_to_include)
{

     using namespace genfile::bgen ;
     using namespace Rcpp;
     m_isQuery = false;
     if(t_bgenFileIndex == ""){
        m_isQuery = false;
        std::cout << "no index file for bgen is provided" << std::endl;
     }else if(t_ranges_to_include.nrow() == 0 &&
        t_ids_to_include.size() == 0 &&
        std::cout << "no query list is provided" << std::endl;
        m_isQuery = false;
     }else{
        m_isQuery = true;
     }  
  

      if(m_isQuery){
        m_genoToTest_bgenDosage = View::create( t_bgenFileName ) ;
        IndexQuery::UniquePtr query = IndexQuery::create( t_bgenFileIndex ) ;
	if (t_ranges_to_include.nrow() > 0){
          StringVector const& chromosome = t_ranges_to_include["chromosome"] ;
          IntegerVector const& start = t_ranges_to_include["start"] ;
          IntegerVector const& end = t_ranges_to_include["end"] ;
          for( int i = 0; i < t_ranges_to_include.nrows(); ++i ) {
            if( end[i] < start[i] ) {
              throw std::invalid_argument( "Range (" + chromosome[i] + ":" + atoi( start[i] ) + "-" + atoi( end[i] ) + ") is malformed." ) ;
             }
             query->include_range( IndexQuery::GenomicRange( std::string( chromosome[i] ), start[i], end[i] )) ;
           }
        }

	if (t_ids_to_include.size() != 0){
          query->include_rsids(t_ids_to_include);
        }

	query->initialise() ;
        if(query->number_of_variants() > 0){
          m_genoToTest_bgenDosage->set_query( query ) ;
          m_M0 = m_genoToTest_bgenDosage->number_of_variants() ;
          m_N0 = m_genoToTest_bgenDosage->number_of_samples();
          std::cout << m_M0 << " markers will be analyzed " << std::endl;
          //return numMarkers ;
        }else{
          std::cout << "No queried variant is found in the bgen file! All variants bgen file will be analyzed" << std::endl;
          m_isQuery = false;
        }
      }else{
          gm_stream.reset(
            new std::ifstream( t_bgenFileName.c_str(), std::ifstream::binary )
          ) ;

          if( !*gm_stream ) {
            throw std::invalid_argument( t_bgenFileName ) ;
          }
	  
	  gm_stream->seekg( 0, std::ios::beg ) ;
          genfile::bgen::read_offset( *gm_stream, &gm_offset ) ;
          //printf("2\n");fflush(NULL);
          genfile::bgen::read_header_block( *gm_stream, &gm_context ) ;

          m_N0 = gm_context.number_of_samples;
          numSamples_bgen = int(m_N0);
          std::cout << numSamples_bgen << " samples are found in the bgen file" << std::endl;

          // Jump to the first variant data block.
          gm_stream->seekg( gm_offset + 4 ) ;
          //printf("4\n");fflush(NULL);
          uint Mbgen = gm_context.number_of_variants;
          m_M0 = int(Mbgen);
          //std::cout << "All " << numMarkers << " markers will be analyzed " << std::endl;
          std::cout << m_M0 << " markers are found in the bgen file " << std::endl;
       }
  }



  /****code from BOLT-LMM v2.3.4***/
  
  /********** READ HEADER v1.2**********/
  //m_ibgenFile.open(m_bgenFile.c_str(), std::ios::binary);
/* 
  m_fin = fopen(t_bgenFileName.c_str(), "rb");
  //chr offset; m_ibgenFile.read(&offset, 4, 1); //cout << "offset: " << offset << endl;
  uint32_t offset; fread(&offset, 4, 1, m_fin); //cout << "offset: " << offset << endl;
  uint32_t L_H; fread(&L_H, 4, 1, m_fin); //cout << "L_H: " << L_H << endl;
  
  
  fread(&m_M0, 4, 1, m_fin); std::cout << "snpBlocks (Mbgen): " << m_M0 << std::endl;
  assert(m_M0 != 0);
  //unsigned int Nbgen; fread(&Nbgen, 4, 1, m_fin); std::cout << "samples (Nbgen): " << Nbgen << std::endl;
  fread(&m_N0, 4, 1, m_fin); std::cout << "samples (Nbgen): " << m_N0 << std::endl;
  unsigned int m_Nsample = t_SampleInBgen.size();
  m_SampleInBgen = t_SampleInBgen;
  if (m_N0 != m_Nsample) {
    std::cerr << "ERROR: Number of samples in BGEN header does not match sample file" << std::endl;
    exit(1);
  }
  char magic[5]; fread(magic, 1, 4, m_fin); magic[4] = '\0'; //cout << "magic bytes: " << string(magic) << endl;
  fseek(m_fin, L_H-20, SEEK_CUR); //cout << "skipping L_H-20 = " << L_H-20 << " bytes (free data area)" << endl;
  uint32_t flags; fread(&flags, 4, 1, m_fin); //cout << "flags: " << flags << endl;
  uint32_t CompressedSNPBlocks = flags&3; std::cout << "CompressedSNPBlocks: " << CompressedSNPBlocks << std::endl;
  assert(CompressedSNPBlocks==1); // REQUIRE CompressedSNPBlocks==1
  uint32_t Layout = (flags>>2)&0xf; std::cout << "Layout: " << Layout << std::endl;
  assert(Layout==1 || Layout==2); // REQUIRE Layout==1 or Layout==2
  fseek(m_fin, offset+4, SEEK_SET);
}
*/

void BgenClass::setPosSampleInBgen(std::vector<std::string> & t_SampleInModel)
{
  std::cout << "Setting position of samples in Bgen files...." << std::endl;	  
  m_N = t_SampleInModel.size();
  
  // updated by BWJ on 03/14/2021
  
  Rcpp::CharacterVector SampleInBgen(m_N0);
  for(uint32_t i = 0; i < m_N0; i++)
    SampleInBgen(i) = m_SampleInBgen.at(i);
  
  Rcpp::CharacterVector SampleInModel(m_N);
  for(uint32_t i = 0; i < m_N; i++)
    SampleInModel(i) = t_SampleInModel.at(i);
  
  Rcpp::IntegerVector posSampleInBgen = Rcpp::match(SampleInModel, SampleInBgen);
  for(uint32_t i = 0; i < m_N; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInBgen.at(i)))
      Rcpp::stop("At least one subject requested is not in BGEN file.");
  }
  
  Rcpp::IntegerVector posSampleInModel = Rcpp::match(SampleInBgen, SampleInModel);
  m_posSampleInModel.resize(m_N0);
  for(uint32_t i = 0; i < m_N0; i++){
    if(Rcpp::IntegerVector::is_na(posSampleInModel.at(i))){
      m_posSampleInModel.at(i) = -1;
    }else{
      m_posSampleInModel.at(i) = posSampleInModel.at(i) - 1;   // convert "starting from 1" to "starting from 0"
    }
    //std::cout << "m_posSampleInModel.at(i) " << i << " " << m_posSampleInModel.at(i) << std::endl; 
  }
  
  // end of the update on 03/14/2021
}


void BgenClass::Parse2(unsigned char *buf, uint bufLen, const unsigned char *zBuf, uint zBufLen,std::string & snpName,arma::vec & dosages, double & AC, double & AF, std::vector<uint> & indexforMissing, double & info, std::vector<uint> & indexNonZero, bool isImputation) {
 
  //arma::vec timeoutput0 = getTime();

  size_t destLen = bufLen;

  unsigned char *bufAt = buf;
  uint N = bufAt[0]|(bufAt[1]<<8)|(bufAt[2]<<16)|(bufAt[3]<<24); bufAt += 4;
  
  if (N != m_N0) {
    std::cerr << "ERROR: " << snpName << " has N = " << N << " (mismatch with header block)" << std::endl;
    exit(1);
  }
  uint K = bufAt[0]|(bufAt[1]<<8); bufAt += 2;
  if (K != 2U) {
    std::cerr << "ERROR: " << snpName << " has K = " << K << " (non-bi-allelic)" << std::endl;
    exit(1);
  }
  uint Pmin = *bufAt; bufAt++;
  if (Pmin != 2U) {
    std::cerr << "ERROR: " << snpName << " has minimum ploidy = " << Pmin << " (not 2)" << std::endl;
    exit(1);
  }
  uint Pmax = *bufAt; bufAt++;
  if (Pmax != 2U) {
    std::cerr << "ERROR: " << snpName << " has maximum ploidy = " << Pmax << " (not 2)" << std::endl;
    exit(1);
  }
  
  
  const unsigned char *ploidyMissBytes = bufAt;
//  std::cout << "N " << N << std::endl;
  for (uint i = 0; i < N; i++) {
    uint ploidyMiss = *bufAt; bufAt++;
    if (ploidyMiss != 2U && ploidyMiss != 130U) {
      std::cerr << "ERROR: " << snpName << " has ploidy/missingness byte = " << ploidyMiss
                << " (not 2 or 130)" << std::endl;
      exit(1);
    }
  }



  uint Phased = *bufAt; bufAt++;
  if (Phased != 0U) {
    std::cerr << "ERROR: " << snpName << " has Phased = " << Pmax << " (not 0)" << std::endl;
    exit(1);
  }
  uint B = *bufAt; bufAt++;
  if (B != 8U) {
    std::cerr << "ERROR: " << snpName << " has B = " << B << " (not 8)" << std::endl;
    exit(1);
  }
  double lut[256];
  for (int i = 0; i <= 255; i++)
    lut[i] = i/255.0;

  //arma::vec timeoutput2 = getTime();

    double sum_eij = 0, sum_fij_minus_eij2 = 0, sum_eij_sub = 0, sum_fij_minus_eij2_sub = 0; // for INFO
    double p11,p10,p00,dosage,eij,fij, eijsub, fijsub;
    dosages.clear();
    //dosages.resize(m_N);
    dosages.set_size(m_N);
    dosages.fill(arma::datum::nan);
    if(!m_isSparseDosagesInBgen){
      //dosages.resize(m_N);
      dosages.set_size(m_N);
    }
    std::size_t missing_cnt = 0;
    double dosage_new;
    for (uint i = 0; i < N; i++) {
     //if(i == 1){std::cout << "ploidyMissBytes[i] " << ploidyMissBytes[i] << std::endl;}
     if (ploidyMissBytes[i] != 130U){
      //bufAt += 2;
      p11 = lut[*bufAt]; bufAt++;
      p10 = lut[*bufAt]; bufAt++;
      uint ninzeroind = 0;
     if(m_posSampleInModel[i] >= 0){
      p00 = 1 - p11 - p10; //can remove
      dosage = 2*p11 + p10;

           dosage_new = 2-dosage;

        eij = dosage;
     //if(m_posSampleInModel[i] >= 0){
        fij = 4*p11 + p10;
        sum_eij += eij;
        sum_fij_minus_eij2 += fij - eij*eij;
          if(!m_isSparseDosagesInBgen){
              dosages[m_posSampleInModel[i]] = dosage_new;
              if(dosage_new > 0){
                indexNonZero.push_back(m_posSampleInModel[i]);
              }
          }else{
              if(dosage_new > 0){
                dosages[i] = dosage_new;
		ninzeroind = ninzeroind + 1;      
                indexNonZero.push_back(m_posSampleInModel[i]);
              }
          }
          sum_eij_sub += eij;
      }



     }else if(ploidyMissBytes[i] == 130U){
        bufAt += 2;
        if(m_posSampleInModel[i] >= 0){
          indexforMissing.push_back(m_posSampleInModel[i]);
          ++missing_cnt;
          if(!m_isSparseDosagesInBgen){
            dosages[m_posSampleInModel[i]] = -1;
          }
        }
     }
    }	  

   //arma::vec timeoutput3 = getTime();
   //printTime(timeoutput2, timeoutput3, "Parse2");
    //}
  //std::cout << "sum_eij_sub: " << sum_eij_sub << std::endl;
  
     //if(m_AlleleOrder == "alt-first"){
     //      AC = sum_eij_sub;
     // }else{
      	   AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;
            //dosage = 2-dosage;
      //}
   
  //AC = 2* ((double) (m_N - missing_cnt)) - sum_eij_sub;
  if(m_N == missing_cnt){
    AF = 0;
  }else{
    AF = AC/ 2/ ((double) (m_N - missing_cnt)) ;
  }
  
  //std::cout << "AC: " << AC << std::endl;
  double thetaHat = sum_eij / (2* (m_N - missing_cnt));
  //std::cout << "missing_cnt " << sum_eij << std::endl;
  if(isImputation){
    info = thetaHat==0 || thetaHat==1 ? 1 :
    1 - sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat));
    //std::cout << "info " << info << std::endl;
    //std::cout << "sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat)); " << sum_fij_minus_eij2 / (2*(m_N - missing_cnt)*thetaHat*(1-thetaHat)) << std::endl;
  }else{
	info=1.0;
  }	  
//arma::vec timeoutput4 = getTime();
//}
//printTime(timeoutput0, timeoutput1, "time 0 to 1 Parse2");
//printTime(timeoutput1, timeoutput2, "time 1 to 2 Parse2");
//printTime(timeoutput2, timeoutput3, "time 2 to 3 Parse2");
//printTime(timeoutput3, timeoutput4, "time 3 to 4 Parse2"); 
//printTime(timeoutput3, timeoutput4, "time 3 to 4 Unified_getOneMarker");

}

void BgenClass::getOneMarker(
                                  std::string& t_ref,   // REF allele
                                  std::string& t_alt,   // ALT allele (should probably be minor allele, otherwise, computation time will increase)
                                  std::string& t_marker,    // marker ID extracted from genotype file
                                  uint32_t& t_pd,           // base position
                                  std::string& t_chr,       // chromosome
                                  double& t_altFreq,        // frequency of ALT allele
                                  double& t_altCounts,      // counts of ALT allele
                                  double& t_missingRate,    // missing rate
                                  double& t_imputeInfo,     // imputation information score, i.e., R2 (all 1 for PLINK)
                                  bool & t_isOutputIndexForMissing,         // if true, output index of missing genotype data
                                  std::vector<uint>& t_indexForMissing,     // index of missing genotype data
                                  bool & t_isOnlyOutputNonZero,             // is true, only output a vector of non-zero genotype. (NOTE: if ALT allele is not minor allele, this might take much computation time)
                                  std::vector<uint>& t_indexForNonZero,
                                  bool & t_isBoolRead,        // only used in BGEN, Wei, if you want, you can add details here.
				  arma::vec & dosages,
				  bool t_isImputation)
{

  using namespace genfile;
  using namespace Rcpp ;

  bool temp;
  std::string SNPID, RSID, chromosome, first_allele,second_allele ;
  genfile::bgen::uint32_t position;
  std::vector< std::string > alleles ;
  std::vector< std::vector< double > > probs ;
  std::vector< byte_t > buffer1;
  std::vector< byte_t > buffer2;
  //std::vector< double > dosages;
  double AC, AF;

  m_isReadVariantBgen = genfile::bgen::read_snp_identifying_data(
                        *gm_stream,
                        gm_context,
                        &SNPID,
                        &RSID,
                        &chromosome,
                        &position,
                        &first_allele,
                        &second_allele
                ) ;
  genfile::bgen::read_genotype_data_block(
                        *gm_stream,
                        gm_context,
                        &buffer1
                ) ;

  genfile::bgen::uncompress_probability_data(
                        gm_context,
                        buffer1,
                        &buffer2
                ) ;

  unsigned char * buf  = (unsigned char *) buffer2.data();
  uint Nbgen = gm_context.number_of_samples;
  
  std::string SNPID, RSID, chromosome, first_allele,second_allele ;
  uint32_t position;
  std::vector< std::string > alleles ;
  double AC, AF, info;
  t_indexForMissing.clear();
  t_indexForNonZero.clear();

   Parse2(buf, buffer2.size(), RSID, dosages, AC, AF, t_indexForMissing, info, t_indexForNonZero, t_isImputation); 
    
    t_alt = second_allele;  // default setting is "ref-first" (03-14-2021)
    t_ref = first_allele;
    t_marker = RSID;    // marker ID extracted from genotype file
    t_pd = position;           // base position
    t_chr = chromosome;       // chromosome
    t_altFreq = AF;        // frequency of ALT allele
    t_altCounts = AC;      // counts of ALT allele
    t_imputeInfo = info;     // imputation information score, i.e., R2 (all 1 for PLINK)
    t_missingRate = (double) t_indexForMissing.size() / (double) m_N;    // missing rate
    //t_indexForNonZero = indexNonZero;   
    //std::cout << "indexNonZero.size() " << indexNonZero.size() << std::endl;
    //std::cout << "t_indexForNonZero.size() " << t_indexForNonZero.size() << std::endl;
    
  //arma::vec getOneMarker_timeoutput3 = getTime();

   //printTime(getOneMarker_timeoutput2,getOneMarker_timeoutput3, "getOneMarker_timeoutput2_3");

   /*
    // dosages.clear();
  // indexforMissing.clear();
  //return dosages;

  */
}


void BgenClass::setIsSparseDosageInBgen (bool t_isSparseDosageInBgen){
  m_isSparseDosagesInBgen = t_isSparseDosageInBgen;
}

void BgenClass::setMarkerIndicesToIncludeInBgen (std::vector< int > & t_markerIndicesToInclude){
  m_markerIndicesToInclude = t_markerIndicesToInclude;
}

void BgenClass::setIsDropMissingDosagesInBgen (bool t_isDropmissingdosagesInBgen){
  m_isDropMissingDosagesInBgen = t_isDropmissingdosagesInBgen;
}

void BgenClass::closegenofile(){
	fclose(m_fin);
}

}  
