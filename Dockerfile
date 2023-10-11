ARG BASE_IMAGE=ubuntu:22.04
FROM $BASE_IMAGE as builder

# install binaries needed to compile jonathonl/shrinkwrap & statgen/savvy
# install superlu-dev binary, avoiding install via cget
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    apt-utils \
    cmake \
    g++ \
    liblzma-dev \
    libsuperlu-dev \
    libzstd-dev \
    python3 \
    python3-venv \
    zlib1g-dev

# in Ubuntu 23.04 and later, an error is now passed for attempting to
# globally install packages with an external package manager like pip.
# to solve: use virtual environments for package installs
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# install dependencies (wheel first to avoid future deprecation) :
RUN pip install wheel 
RUN pip install click six cget

# copy over SAIGE repository & remove current install configuration
# instead, cget install below to utilize ubuntu binaries & avoid duplication
WORKDIR /app
COPY . .
RUN rm configure
RUN rm -rf thirdParty

# ignore dependencies from jonathon1/shrinkwrap :
# xz, zlib, zstd are picked up above via binaries
# xz : liblzma-dev, zlib : zlib1g-dev, zstd : libzstd-dev
RUN cget ignore --prefix thirdParty/cget xz zlib zstd
RUN cget install --prefix thirdParty/cget jonathonl/shrinkwrap
RUN cget install --prefix thirdParty/cget statgen/savvy

# install all possible R dependencies from pre-packaged ubuntu binaries
# - much faster install than 'install.packages()' inside R
RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    r-base \
    r-cran-biocmanager \
    r-bioc-biocversion \
    r-cran-devtools \
    r-cran-optparse \
    r-cran-qlcmatrix \
    r-cran-r.utils \
    r-cran-rcppparallel \
    r-cran-roxygen2 \
    r-cran-rversions

# pick up remaining bioconductor dependencies not packaged w/ ubuntu 
RUN Rscript -e 'BiocManager::install("SPAtest")'
RUN Rscript -e 'BiocManager::install("RhpcBLASctl")'
RUN Rscript -e 'BiocManager::install("SKAT")'
RUN Rscript -e 'BiocManager::install("MetaSKAT")'

# install SAIGE, will default into /usr/local/lib/R
RUN R CMD INSTALL .

# ------

# with SAIGE installed, rebuild container with only the runtime dependencies to
# achieve a more lightweight final build
FROM $BASE_IMAGE

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    r-base \
    r-cran-optparse \
    r-cran-rcppparallel

# copy SAIGE run scripts from builder & put into /usr/local/bin for persistence
COPY --from=builder /app/extdata/*.R /usr/local/bin/
RUN chmod ug+x /usr/local/bin/*.R
# copy compiled site-packages from builder
COPY --from=builder /usr/local/lib/R /usr/local/lib/R

RUN createSparseGRM.R  --help
RUN step1_fitNULLGLMM.R --help
RUN step2_SPAtests.R --help
RUN step3_LDmat.R --help

