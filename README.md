HERE we are maintaining an newly improved stable version of SAIGE and SAIGE-GENE+. 
Please find the https://saigegit.github.io/SAIGE-doc/ for documentation.


SAIGE is an R package developed with Rcpp for genome-wide association tests in large-scale data sets and biobanks. The method

- accounts for sample relatedness based on the generalized mixed models
- allows for model fitting with either full or sparse genetic relationship matrix (GRM)
- works for quantitative and binary traits
- handles case-control imbalance of binary traits
- computationally efficient for large data sets
- performs single-variant association tests
- provides effect size estimation through Firth's Bias-Reduced Logistic Regression
- performs conditional association analysis

SAIGE-GENE (now known as SAIGE-GENE+) are new method extension in the R package for testing rare variant in set-based tests.
- performs BURDEN, SKAT, and SKAT-O tests
- allows for tests on multiple minor allele frequencies cutoffs and functional annotations
- allows for specifying weights for markers in the set-based tests
- performs conditional analysis to identify associations independent from nearly GWAS signals


The package takes genotype file input in the following formats
- PLINK (bed, bim, fam), PGEN, BGEN, VCF, BCF, SAV

How to install:

On Linux:
```
# install pixi
curl -fsSL https://pixi.sh/install.sh | sh
# --- need a new shell or source .*shrc here ---

# install deps from pixi.toml
pixi install

# install extra dep
pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'

# set up plink2
curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx
mv plink-ng-2.0.0-a.6.16 plink-ng
pixi run x86_64-conda-linux-gnu-cc -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o plink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate
mv plink2_includes.a .pixi/envs/default/lib

# install SAIGE
pixi run R CMD INSTALL .

# example run
pixi run Rscript extdata/step1_fitNULLGLMM.R --help
```

On MacOS
```
# install pixi
curl -fsSL https://pixi.sh/install.sh | sh
# --- need a new shell or source .*shrc here ---

# install deps from pixi.toml
pixi install

# install extra dep
pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'

# set up plink2
curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx
mv plink-ng-2.0.0-a.6.16 plink-ng
## pixi run clang++ -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o plink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate
# mv plink2_includes.a .pixi/envs/default/lib
mkdir -p plink-ng/2.0/build
cd plink-ng/2.0
### compile using clang
pixi run clang++ -std=c++14 -fPIC -O3 -I../../../.pixi/envs/default/include -c include/*.cc
mv *.o build/
cd build
ar rcs ../plink2_includes.a *.o

cp ../plink2_includes.a ../../../.pixi/envs/default/lib/libplink2_includes.a
cd ../../../
# install SAIGE
pixi run R CMD INSTALL .

# example run
pixi run Rscript extdata/step1_fitNULLGLMM.R --help
```


