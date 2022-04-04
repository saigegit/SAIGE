void writeOutfile_BURDEN(std::string regionName,
			std::vector<std::string>  & BURDEN_AnnoName_Vec,
			std::vector<std::string> & BURDEN_maxMAFName_Vec,
			arma::vec & BURDEN_pval_Vec,
			std::vector<double> & BURDEN_Beta_Vec,
			std::vector<double> & BURDEN_seBeta_Vec,
			arma::vec & BURDEN_pval_cVec,
			std::vector<double> & BURDEN_Beta_cVec,
			std::vector<double> & BURDEN_seBeta_cVec,
			arma::vec & MAC_GroupVec,
			arma::vec & MACCase_GroupVec,
			arma::vec & MACControl_GroupVec,
			arma::vec & NumRare_GroupVec,
			arma::vec & NumUltraRare_GroupVec,
			double cctpval,
			double cctpval_cond,
			unsigned int q_anno,
			unsigned int q_maf,
			bool isCondition,
			std::string t_traitType){

     for(unsigned int j = 0; j < q_anno; j++){
       for(unsigned int m = 0; m < q_maf; m++){
           jm = j*q_maf+m;
           i = jm;	
           OutFile << regionName;
           OutFile << "\t";
           OutFile << BURDEN_AnnoName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_maxMAFName_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_pval_Vec(i);
           OutFile << "\t";
           OutFile << BURDEN_Beta_Vec.at(i);
           OutFile << "\t";
           OutFile << BURDEN_seBeta_Vec.at(i);
           OutFile << "\t";
           if(isCondition){
               OutFile << BURDEN_pval_cVec(i);
               OutFile << "\t";
               OutFile << BURDEN_Beta_cVec.at(i);
               OutFile << "\t";
               OutFile << BURDEN_seBeta_cVec.at(i);
               OutFile << "\t";
	   }
	   OutFile << MAC_GroupVec(i);
           OutFile << "\t";
           if(t_traitType == "binary"){
               OutFile << MACCase_GroupVec(i);
               OutFile << "\t";
               OutFile << MACControl_GroupVec(i);
               OutFile << "\t";
           }
           OutFile << NumRare_GroupVec(i);
           OutFile << "\t";
           OutFile << NumUltraRare_GroupVec(i);
           OutFile << "\n";
	}
     }
     OutFile << regionName;
     OutFile << "\tCauchy\tNA\t";
     OutFile << cctpval;
     OutFile << "\tNA\tNA\t";	
     if(isCondition){
	OutFile << cctpval_cond;
	OutFile << "\tNA\tNA\t";
     }
     OutFile << "NA\t";
     if(t_traitType == "binary"){
        OutFile << "NA\t";
        OutFile << "NA\t";
     }
     OutFile << "NA\t";
     OutFile << "NA\n";	
}
