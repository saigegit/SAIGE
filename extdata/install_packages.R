#!/usr/bin/env Rscript
#install required R packages, from Finnge/SAIGE-IT

#req_packages <- c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "MetaSKAT", "roxygen2", "rversions","devtools", "SKAT")
req_packages <- c("R.utils", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "roxygen2", "rversions","devtools", "SKAT", "RhpcBLASctl", "qlcMatrix", "dplyr")
for (pack in req_packages) {
    if(!require(pack, character.only = TRUE)) {
        install.packages(pack, repos = "https://cloud.r-project.org")
        print(packageVersion(pack))
    }
}
devtools::install_version("Rcpp", version = "1.0.7", repos = "http://cran.us.r-project.org")
print(packageVersion("Rcpp"))
install.packages("RSQLite", repos='http://cran.rstudio.com/')
print(packageVersion("RSQLite"))


#devtools::install_github("leeshawn/SPAtest")
#devtools::install_github("leeshawn/MetaSKAT")
#devtools::install_github("leeshawn/SKAT")
