#!/usr/bin/env Rscript
#install required R packages, from Finnge/SAIGE-IT

install.packages("ellipsis", repos='http://cran.rstudio.com/')
print(packageVersion("ellipsis"))
install.packages("vctrs", repos='http://cran.rstudio.com/')
print(packageVersion("vctrs"))


req_packages <- c("R.utils", "devtools", "Rcpp", "RcppParallel", "RcppArmadillo", "data.table", "RcppEigen", "Matrix", "methods", "BH", "optparse", "SPAtest", "roxygen2", "rversions","devtools", "SKAT", "RhpcBLASctl", "qlcMatrix", "dplyr", "dbplyr")
for (pack in req_packages) {
    if(!require(pack, character.only = TRUE)) {
        install.packages(pack, repos = "https://cloud.r-project.org", dependencies=TRUE)
        print(packageVersion(pack))
    }
}

github_packages <- c("leeshawn/MetaSKAT")
for (pack in github_packages) {
    if(!require(pack, character.only = TRUE)) {
        devtools::install_github(pack)
    }
}

if((!require("RcppArmadillo", character.only = TRUE)) | packageVersion("RcppArmadillo") < "0.10.7.5.0"){
    install.packages("RcppArmadillo", repos = "https://cloud.r-project.org", dependencies=TRUE)
}


devtools::install_version("Rcpp", version = "1.0.7", repos = "http://cran.us.r-project.org")
print(packageVersion("Rcpp"))
install.packages("RSQLite", repos='http://cran.rstudio.com/')
print(packageVersion("RSQLite"))

