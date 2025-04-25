FROM ubuntu:20.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ENV R_VERSION=4.2.3 \
    DEBIAN_FRONTEND=noninteractive \
    PKG_CONFIG_PATH=/usr/lib/pkgconfig:/usr/share/pkgconfig


RUN apt-get update -qq && \
    apt-get -y install --no-install-recommends \
    ca-certificates \
    build-essential \
    gfortran \
    libreadline-dev \
    xorg-dev \
    libbz2-dev \
    liblzma-dev \
    curl \
    libxml2-dev \
    libcairo2-dev \
    libsqlite3-dev \
    libmariadbd-dev \
    libpq-dev \
    libssh2-1-dev \
    unixodbc-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libsodium-dev \
    libtool \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libboost-all-dev \
    pkg-config \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Download and install R
RUN wget -c https://cran.r-project.org/src/base/R-4/R-${R_VERSION}.tar.gz \
    && tar -xf R-${R_VERSION}.tar.gz \
    && cd R-${R_VERSION} \
    && ./configure \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf R-${R_VERSION} R-${R_VERSION}.tar.gz

RUN echo "here"

# Set the CRAN repository and install R.utils and devtools
#RUN R -e "options(repos = c(CRAN = 'https://cloud.r-project.org')); install.packages(c('R.utils', 'devtools'))"

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive apt-get install -y \
    build-essential \
    cmake \
    libopenblas-base \
    python3-pip

RUN pip3 install cget

WORKDIR /app
COPY . .

#RUN Rscript extdata/install_packages_devtools.R

RUN Rscript extdata/install_packages.R

# Force step_2 to use 1 single thread. More threads are ineffective
ENV OMP_NUM_THREADS=1

RUN apt-get update && apt-get install -y \
    liblapack-dev \
    && apt-get clean

# Set environment variables for LAPACK
ENV LD_LIBRARY_PATH="/usr/lib:${LD_LIBRARY_PATH}"
ENV PKG_CONFIG_PATH="/usr/lib/pkgconfig:/usr/share/pkgconfig:${PKG_CONFIG_PATH}"



RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/BH/BH_1.78.0-0.tar.gz', repos=NULL, type='source')"


RUN R CMD INSTALL .
RUN echo "here"
RUN mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/

RUN chmod a+x /usr/local/bin/step1_fitNULLGLMM.R
RUN chmod a+x /usr/local/bin/step2_SPAtests.R
RUN chmod a+x /usr/local/bin/step3_LDmat.R
RUN chmod a+x /usr/local/bin/createSparseGRM.R

RUN createSparseGRM.R  --help
RUN step1_fitNULLGLMM.R --help
RUN step2_SPAtests.R --help
RUN step3_LDmat.R --help

RUN apt-get update
RUN apt-get install time
