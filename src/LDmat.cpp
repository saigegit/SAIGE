// This includes the main codes to connect C++ and R

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include <vector>
#include <thread>         // std::this_thread::sleep_for
#include <chrono>         // std::chrono::seconds
// std::this_thread::sleep_for (std::chrono::seconds(1));
#include <cstdio>         // std::remove
#include <fstream>
#include <string.h>
// Currently, omp does not work well, will check it later
// error: SET_VECTOR_ELT() can only be applied to a 'list', not a 'character'
// remove all Rcpp::List to check if it works
// #include <omp.h>
// // [[Rcpp::plugins(openmp)]]]

#include "LDmat.hpp"
#include "Main.hpp"
#include "PLINK.hpp"
#include "BGEN.hpp"
#include "VCF.hpp"
#include "SAIGE_test.hpp"
#include "UTIL.hpp"
#include "CCT.hpp"

#include <Rcpp.h>
#include "getMem.hpp"

#include <boost/math/distributions/beta.hpp>

// global objects for different genotype formats

static PLINK::PlinkClass* ptr_gPLINKobj = NULL;
static BGEN::BgenClass* ptr_gBGENobj = NULL;
static VCF::VcfClass* ptr_gVCFobj = NULL;


std::ofstream OutFile_single_LDmat;
std::string g_outputFilePrefixSingle_LDmat;
std::ofstream OutFile_LDmat;
std::string g_outputFilePrefix_LDmat;

// global variables for analysis
std::string g_impute_method_LDmat;      // "mean", "minor", or "drop", //drop is not allowed
double g_dosage_zerod_MAC_cutoff_LDmat;
double g_dosage_zerod_cutoff_LDmat;

double g_missingRate_cutoff_LDmat;
double g_maxMAFLimit_LDmat;
//unsigned int g_omp_num_threads;
double g_marker_minMAF_cutoff_LDmat;
double g_marker_minMAC_cutoff_LDmat;
double g_marker_minINFO_cutoff_LDmat;

unsigned int g_region_maxMarkers_cutoff_LDmat;   // maximal number of markers in one chunk, only used for region-based analysis to reduce memory usage


// [[Rcpp::export]]
void setGlobalVarsInCPP_LDmat(std::string t_impute_method,
                               double t_dosage_zerod_cutoff,
                               double t_dosage_zerod_MAC_cutoff,

                	       double t_missing_cutoff,
			       double t_maxMAFLimit,

                               double t_min_maf_marker,
                               double t_min_mac_marker,
                               double t_min_info_marker,
			       unsigned int t_max_markers_region, 
			       
			       std::string t_outputFile)
{
  g_impute_method_LDmat = t_impute_method;
  g_dosage_zerod_cutoff_LDmat = t_dosage_zerod_cutoff;
  g_dosage_zerod_MAC_cutoff_LDmat = t_dosage_zerod_MAC_cutoff;

  g_missingRate_cutoff_LDmat = t_missing_cutoff;
  g_maxMAFLimit_LDmat = t_maxMAFLimit;

  g_marker_minMAF_cutoff_LDmat = t_min_maf_marker;
  g_marker_minMAC_cutoff_LDmat = t_min_mac_marker;
  g_marker_minINFO_cutoff_LDmat = t_min_info_marker;
  g_region_maxMarkers_cutoff_LDmat = t_max_markers_region;

  g_outputFilePrefixSingle_LDmat = t_outputFile + ".marker_info.txt";

}


// [[Rcpp::export]]
void LDmatRegionInCPP(
                           std::string t_genoType,     // "PLINK", "BGEN"
                           std::vector<std::string> & t_genoIndex_prev,
                           std::vector<std::string> & t_genoIndex,
                           arma::mat & annoIndicatorMat,
                           std::string t_outputFile,
			   unsigned int t_n,           // sample size
                           bool t_isImputation,
                           std::vector<std::string> & annoStringVec,
                           std::string regionName
			   ) 
{
  //create the output list
  //arma::vec timeoutput1 = getTime();
  unsigned int q0 = t_genoIndex.size();                 // number of markers (before QC) in one region
  unsigned int q = q0;
  std::vector<std::string> markerVec(q);
  std::vector<std::string> chrVec(q);  // marker IDs
  std::vector<std::string> posVec(q);  // marker IDs
  std::vector<std::string> refVec(q);  // marker IDs
  std::vector<std::string> altVec(q);  // marker IDs
  std::vector<double> altCountsVec(q, arma::datum::nan);    // allele counts of ALT allele.
  std::vector<double> imputationInfoVec(q, arma::datum::nan);    // imputation info of ALT allele.
  std::vector<double> missingRateVec(q, arma::datum::nan);
  std::vector<double> altFreqVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<std::string> infoVec(q);    // marker information: CHR:POS:REF:ALT
  std::vector<double> MACVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<double> MAFVec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.
  std::vector<uint32_t> N_Vec(q, arma::datum::nan);      // allele frequencies of ALT allele, this is not always < 0.5.

  unsigned int m1 = g_region_maxMarkers_cutoff_LDmat;     // number of markers in all chunks expect for the last chunk
  std::vector<unsigned int> mPassCVVec;
  arma::vec GVec(t_n);
  arma::vec GZeroVec(t_n);
  std::vector<uint> indexZeroVec;
  std::vector<uint> indexNonZeroVec;
  unsigned int nchunks = 0; //number of chunks
  unsigned int ichunk = 0; //ith chunk
  unsigned int i1InChunk = 0; //i1th marker in ith chunk
  unsigned int i1 = 0;    // index of Markers (non-URV)



  unsigned int jm;
  std::vector<uint32_t> location_m_P1Mat;
  std::vector<uint32_t> location_n_P1Mat;
  std::vector<double> value_P1Mat;

  for(unsigned int i = 0; i < q0; i++)
  {
    // marker-level information
    double altFreq, altCounts, missingRate, imputeInfo;
    std::vector<uint32_t> indexForMissing;
    std::string chr, ref, alt, marker;
    uint32_t pd;
    bool flip = false;
    bool isOutputIndexForMissing = true;
    bool isOnlyOutputNonZero = false;
    GVec.resize(t_n);
    GVec.zeros();
    std::string t_genoIndex_str = t_genoIndex.at(i);
    char* end;
    uint64_t gIndex = std::strtoull( t_genoIndex_str.c_str(), &end,10 );
    std::remove(end);
    uint64_t gIndex_prev = 0;

    if(i == 0){
        gIndex_prev = 0;
    }else{
        char* end_prev;
        std::string t_genoIndex_prev_str;
        if(t_genoType == "bgen"){
            t_genoIndex_prev_str = t_genoIndex_prev.at(i-1);
        }else if(t_genoType == "plink"){
            t_genoIndex_prev_str = t_genoIndex.at(i-1);
        }
        gIndex_prev = std::strtoull( t_genoIndex_prev_str.c_str(), &end_prev,10 );
        std::remove(end_prev);
    }

    bool isReadMarker = Unified_getOneMarker(t_genoType, gIndex_prev, gIndex, ref, alt, marker, pd, chr, altFreq, altCounts, missingRate, imputeInfo,
                                          isOutputIndexForMissing, // bool t_isOutputIndexForMissing,
                                          indexForMissing,
                                          isOnlyOutputNonZero, // bool t_isOnlyOutputNonZero,
                                          indexNonZeroVec,
                                          GVec,
                                          t_isImputation);



    if(!isReadMarker){
      std::cout << "ERROR: Reading " <<  i << "th marker failed." << std::endl;
      break;
    }
    std::string pds = std::to_string(pd);
    std::string info = chr+":"+std::to_string(pd)+":"+ref+":"+alt;

    double MAF = std::min(altFreq, 1 - altFreq);
    double w0;
    double MAC = MAF * 2 * t_n * (1 - missingRate);   // checked on 08-10-2021
    flip = imputeGenoAndFlip(GVec, altFreq, altCounts, indexForMissing, g_impute_method_LDmat, g_dosage_zerod_cutoff_LDmat, g_dosage_zerod_MAC_cutoff_LDmat, MAC, indexZeroVec, indexNonZeroVec);
    arma::uvec indexZeroVec_arma, indexNonZeroVec_arma;
    MAF = std::min(altFreq, 1 - altFreq);
    MAC = std::min(altCounts, t_n *2 - altCounts);
    chrVec.at(i) = chr;
    posVec.at(i) = pds;
    if(flip){	
	refVec.at(i) = alt;
	altVec.at(i) = ref;
    }else{	    
        refVec.at(i) = ref;
        altVec.at(i) = alt;
    }
    
    markerVec.at(i) = marker;             // marker IDs
    altFreqVec.at(i) = altFreq;           // allele frequencies of ALT allele, this is not always < 0.5.
    missingRateVec.at(i) = missingRate;
    MAFVec.at(i) = MAF;
    N_Vec.at(i) = t_n;
    imputationInfoVec.at(i) = imputeInfo;
    infoVec.at(i) = info;                 // marker information: CHR:POS:REF:ALT

   if((missingRate > g_missingRate_cutoff_LDmat) || (MAF > g_maxMAFLimit_LDmat) || (MAF < g_marker_minMAF_cutoff_LDmat) || (MAC < g_marker_minMAC_cutoff_LDmat) || (imputeInfo < g_marker_minINFO_cutoff_LDmat)){
	       continue;
   }else{
    MACVec.at(i) = MAC;
    altCountsVec.at(i) = altCounts;
    indexNonZeroVec_arma = arma::conv_to<arma::uvec>::from(indexNonZeroVec);
    location_m_P1Mat.insert(std::end(location_m_P1Mat), indexNonZeroVec.size(), i1InChunk);
    location_n_P1Mat.insert(std::end(location_n_P1Mat), std::begin(indexNonZeroVec), std::end(indexNonZeroVec));
    indexNonZeroVec.clear();
    arma::vec GVecnonZero = GVec(indexNonZeroVec_arma);
    typedef std::vector<double> stdvec;
    stdvec GVecnonZero_stdvec = arma::conv_to< stdvec >::from(GVecnonZero);
    value_P1Mat.insert(std::end(value_P1Mat), std::begin(GVecnonZero_stdvec), std::end(GVecnonZero_stdvec));   
    GVecnonZero_stdvec.clear(); 
    i1 += 1;
    i1InChunk += 1;
   }//continue; else


     if(i1InChunk == m1 ){
        std::cout << "In chunks 0-" << ichunk << ", " <<  i1 << " markers are included." << std::endl;
	int nonzeroSize = location_m_P1Mat.size();
	arma::uvec location_m_P1Mat_arma =  arma::conv_to<arma::uvec>::from(location_m_P1Mat);
	arma::uvec location_n_P1Mat_arma =  arma::conv_to<arma::uvec>::from(location_n_P1Mat);

	arma::umat location_P1Mat_arma (2, location_m_P1Mat_arma.n_elem);
        location_P1Mat_arma.row(0) = location_m_P1Mat_arma.t();
	location_P1Mat_arma.row(1) = location_n_P1Mat_arma.t();
	arma::ivec value_P1Mat_arma =  arma::conv_to<arma::ivec>::from(value_P1Mat);
	arma::sp_imat P1Mat(location_P1Mat_arma, value_P1Mat_arma, i1InChunk, t_n);

	P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
	//P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin.txt", arma::coord_ascii);
	value_P1Mat.clear();
	location_m_P1Mat.clear();
	location_n_P1Mat.clear();
        //P2Mat.save(t_outputFile + "_P2Mat_Chunk_" + std::to_string(ichunk) + ".bin");

        mPassCVVec.push_back(m1);
        ichunk += 1;
        i1InChunk = 0;
        nchunks = nchunks + 1;
    }
    Rcpp::checkUserInterrupt();
  }// for(unsigned int i = 0; i < q0; i++)


  //the second last chunk
  if(i1InChunk != 0){
      std::cout << "In chunks 0-" << ichunk << ", " <<  i1 << " markers are included." << std::endl;
	int nonzeroSize = location_m_P1Mat.size();
        arma::uvec location_m_P1Mat_arma =  arma::conv_to<arma::uvec>::from(location_m_P1Mat);
        arma::uvec location_n_P1Mat_arma =  arma::conv_to<arma::uvec>::from(location_n_P1Mat);
        arma::umat location_P1Mat_arma (2, location_m_P1Mat_arma.n_elem);
        location_P1Mat_arma.row(0) = location_m_P1Mat_arma.t();
        location_P1Mat_arma.row(1) = location_n_P1Mat_arma.t();
       
       	arma::ivec value_P1Mat_arma =  arma::conv_to<arma::ivec>::from(value_P1Mat);
        arma::sp_imat P1Mat(location_P1Mat_arma, value_P1Mat_arma, i1InChunk, t_n);
        P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin");
        //P1Mat.save(t_outputFile + "_P1Mat_Chunk_" + std::to_string(ichunk) + ".bin.txt", arma::coord_ascii);
        value_P1Mat.clear();
        location_m_P1Mat.clear();
        location_n_P1Mat.clear();

    ichunk = ichunk + 1;
    //}
    mPassCVVec.push_back(i1InChunk);
    nchunks = nchunks + 1;
    i1InChunk = 0;
  }


  int mPassCVVecsize = mPassCVVec.size();
  nchunks = mPassCVVecsize;


  //arma::imat VarMat;
  //VarMat.resize(i1, i1);
 
  //arma::uvec rowIndices_VarMat;
  //arma::uvec colIndices_VarMat;
  //arma::ivec values_VarMat;
  // Get the nonzero locations of the first row
  std::vector<unsigned int> rowIndices_VarMat_vec;
  std::vector<unsigned int> colIndices_VarMat_vec;
  unsigned int rowind_new;
  unsigned int colind_new;
  std::vector<int> values_VarMat_vec;

  arma::sp_imat P1Mat; 
  arma::sp_imat P2Mat; 
  // the region includes more markers than limitation, so P1Mat and P2Mat have been put in hard drive
  if(nchunks >= 1)
  {
    int first_row = 0, first_col = 0, last_row = 0, last_col = 0;

    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      last_row = first_row + mPassCVVec.at(index1) - 1;

      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
 
      P1Mat.load(P1MatFile);

       //P1Mat.print("P1Mat");

      if(P1Mat.n_cols == 0) continue;

      // off-diagonal sub-matrix
      for(unsigned int index2 = 0; index2 < index1; index2++)
      {
        std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index2 << "/" << nchunks - 1 << ")........" << std::endl;

        P2Mat.load(t_outputFile + "_P1Mat_Chunk_" + std::to_string(index2) + ".bin");
       //P2Mat.print("P2Mat");

        if(P2Mat.n_cols == 0) continue;
        arma::sp_imat offVarMat = P1Mat * (P2Mat.t());
         
	last_col = first_col + mPassCVVec.at(index2) - 1;

        //VarMat.submat(first_row, first_col, last_row, last_col) = offVarMat;
        //VarMat.submat(first_col, first_row, last_col, last_row) = offVarMat.t();

  //arma::uword nnz = offVarMat.n_nonzero;
  // Create vectors to store row indices, column indices, and values
  //arma::uvec rowIndices(nnz);
  //arma::uvec colIndices(nnz);
  //arma::vec values(nnz);
  // Iterate over non-zero elements and extract information
  arma::sp_imat::const_iterator it = offVarMat.begin();
  arma::sp_imat::const_iterator it_end = offVarMat.end();

  for (arma::uword iof = 0; it != it_end; ++it, ++iof) {
    rowind_new = it.row() + first_row;
    colind_new = it.col() + first_col;
    if(rowind_new >= colind_new){ 
      rowIndices_VarMat_vec.push_back(it.row() + first_row);
      colIndices_VarMat_vec.push_back(it.col() + first_col);
      values_VarMat_vec.push_back(*it);
    }
    //rowIndices(iof) = it.row();
    //colIndices(iof) = it.col();
    //values(iof) = *it;
  }

        first_col = last_col + 1;
 }

      // diagonal sub-matrix
      last_col = first_col + mPassCVVec.at(index1) - 1;
      std::cout << "Analyzing chunks (" << index1 << "/" << nchunks - 1 << ", " << index1 << "/" << nchunks - 1 << ")........" << std::endl;
      //std::cout << "P2Mat.n_cols " << P2Mat.n_cols << std::endl;
      P2Mat.load(t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin");

      arma::sp_imat diagVarMat = P1Mat * (P2Mat.t());

      arma::sp_imat::const_iterator it = diagVarMat.begin();
      arma::sp_imat::const_iterator it_end = diagVarMat.end();

  for (arma::uword iof = 0; it != it_end; ++it, ++iof) {
    rowind_new = it.row() + first_row;
    colind_new = it.col() + first_col;
    if(rowind_new >= colind_new){
      rowIndices_VarMat_vec.push_back(it.row() + first_row);
      colIndices_VarMat_vec.push_back(it.col() + first_col);
      values_VarMat_vec.push_back(*it);
    }
    //rowIndices(iof) = it.row();
    //colIndices(iof) = it.col();
    //values(iof) = *it;
  }

      //VarMat.submat(first_row, first_col, last_row, last_col) = diagVarMat;

      first_row = last_row + 1;
      first_col = 0;
      Rcpp::checkUserInterrupt();
    }

    for(unsigned int index1 = 0; index1 < nchunks; index1++)
    {
      std::string P1MatFile = t_outputFile + "_P1Mat_Chunk_" + std::to_string(index1) + ".bin";
      //std::string P2MatFile = t_outputFile + "_P2Mat_Chunk_" + std::to_string(index1) + ".bin";
      const char* File1 = P1MatFile.c_str();
      //const char* File2 = P2MatFile.c_str();
      std::remove(File1);
      //std::remove(File2);
    }

  }

/*
//write the vectors to a binary file
  std::ofstream file(t_outputFile + "_"+regionName+".txt", std::ios::binary);
  if (!file.is_open()) {
    Rcpp::stop("Error opening file for writing.");
  }
  // Write the sizes of the vectors to the file
  size_t size1 = rowIndices_VarMat_vec.size();
  size_t size2 = colIndices_VarMat_vec.size();
  size_t size3 = values_VarMat_vec.size();

  file.write(reinterpret_cast<const char*>(&size1), sizeof(size_t));
  file.write(reinterpret_cast<const char*>(&size2), sizeof(size_t));
  file.write(reinterpret_cast<const char*>(&size3), sizeof(size_t));

  // Write the vector elements to the file
  file.write(reinterpret_cast<const char*>(rowIndices_VarMat_vec.data()), size1 * sizeof(int));
  file.write(reinterpret_cast<const char*>(colIndices_VarMat_vec.data()), size2 * sizeof(int));
  file.write(reinterpret_cast<const char*>(values_VarMat_vec.data()), size3 * sizeof(int));
*/

g_outputFilePrefix_LDmat = t_outputFile + "_"+regionName+".txt";
OutFile_LDmat.open(g_outputFilePrefix_LDmat.c_str());

  //std::ofstream file(t_outputFile + "_"+regionName+".txt");
  if (OutFile_LDmat.is_open()) {
  for(unsigned int k = 0; k < rowIndices_VarMat_vec.size(); k++){
      OutFile_LDmat << rowIndices_VarMat_vec.at(k);
      OutFile_LDmat << " ";
      OutFile_LDmat << colIndices_VarMat_vec.at(k);
      OutFile_LDmat << " ";
      OutFile_LDmat << values_VarMat_vec.at(k);
      OutFile_LDmat << "\n";
  }
 
 // Close the file
  OutFile_LDmat.close();
}else{
  std::cout << g_outputFilePrefix_LDmat << " is not opened" << std::endl;
}
/*
        arma::uvec location_row_VarMat_arma =  arma::conv_to<arma::uvec>::from(rowIndices_VarMat_vec);
        arma::uvec location_col_VarMat_arma =  arma::conv_to<arma::uvec>::from(colIndices_VarMat_vec);
	arma::umat location_VarMat_arma (2, location_row_VarMat_arma.n_elem);
	location_VarMat_arma.row(0) = location_row_VarMat_arma.t();
	location_VarMat_arma.row(1) = location_col_VarMat_arma.t();
	arma::ivec value_VarMat_arma = arma::conv_to<arma::ivec>::from(values_VarMat_vec);
	arma::sp_imat VarMat(location_VarMat_arma, value_VarMat_arma, i1, i1);

  //VarMat.print("VarMat");
  VarMat.save(t_outputFile + "_"+regionName+".txt", arma::coord_ascii);
  VarMat.clear();
*/
  writeOutfile_single_LDmat(chrVec,
		  posVec,
		  refVec,
		  altVec,
		  infoVec,
		  MACVec,
		  missingRateVec,
		  N_Vec,
		  regionName);
  
}


void writeOutfile_single_LDmat(
                        std::vector<std::string> & chrVec,
                        std::vector<std::string> & posVec,
                        std::vector<std::string> & refVec,
                        std::vector<std::string> & altVec,
			std::vector<std::string> & infoVec,
                        std::vector<double> & altCountsVec,
                        std::vector<double> & missingRateVec,
                        std::vector<uint32_t> & N_Vec,
			std::string regionName
){
	
  int numtest = 0;
  for(unsigned int k = 0; k < chrVec.size(); k++){
    if(!std::isnan(altCountsVec.at(k))){
      OutFile_single_LDmat << chrVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << posVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << refVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << altVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << altCountsVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << N_Vec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << missingRateVec.at(k);
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << regionName;
      OutFile_single_LDmat << "\t";
      OutFile_single_LDmat << numtest;
      OutFile_single_LDmat << "\n";
      numtest = numtest + 1;
    }	    
  }
}


// [[Rcpp::export]]
bool openOutfile_single_LDmat(bool isappend){
      bool isopen;
      if(!isappend){
        OutFile_single_LDmat.open(g_outputFilePrefixSingle_LDmat.c_str());
	isopen = OutFile_single_LDmat.is_open();
	if(isopen){	    
          OutFile_single_LDmat << "CHR\tPOS\tMajor_Allele\tMinor_Allele\tMAC\tN\tMissing_rate\tSet\tIndex\n";
         }
       }else{
         OutFile_single_LDmat.open(g_outputFilePrefixSingle_LDmat.c_str(), std::ofstream::out | std::ofstream::app);
         isopen = OutFile_single_LDmat.is_open();
       }
       return(isopen);
}


// [[Rcpp::export]]
void closeOutfile_single_LDmat(){
  OutFile_single_LDmat.close();
}  
