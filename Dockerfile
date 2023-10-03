FROM ubuntu:23.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    apt-utils \
    build-essential \
    cmake \
    git \
    python3 \
    python3-venv

# in Ubuntu 23.04 and later, an error is now passed for attempting to
# globally install packages with an external package manager like pip.
# to solve: use virtual environments for package installs
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"
# install dependencies:
RUN pip install cget

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

WORKDIR /app
COPY . .

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

RUN R CMD INSTALL .

# move run scripts from SAIGE repo to /usr/local/bin to be picked up
# by $PATH & run from anywhere in the container
RUN mv extdata/step1_fitNULLGLMM.R \
    extdata/step2_SPAtests.R \
    extdata/step3_LDmat.R \
    extdata/createSparseGRM.R \
    /usr/local/bin/

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
RUN chmod a+x /usr/local/bin/step2_SPAtests.R
RUN chmod a+x /usr/local/bin/step3_LDmat.R
RUN chmod a+x /usr/local/bin/createSparseGRM.R

RUN createSparseGRM.R  --help
RUN step1_fitNULLGLMM.R --help
RUN step2_SPAtests.R --help
RUN step3_LDmat.R --help

