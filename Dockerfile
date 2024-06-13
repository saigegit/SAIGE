ARG BASE_IMAGE=ubuntu:24.04
FROM $BASE_IMAGE as base
# ------------------------------------------------------------------------------

# shared builder and runtime dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive \
	apt-get install -y --no-install-recommends --no-install-suggests \
	libsuperlu6 \
    r-base r-base-dev \
	r-cran-data.table r-cran-lattice r-cran-matrix \
	r-cran-optparse r-cran-rcpp r-cran-rcppparallel \
	r-cran-rhpcblasctl \
	&& rm -rf /var/lib/apt/lists/* \
	&& rm -rf /tmp/*
# ------------------------------------------------------------------------------
# BUILDER LAYER
FROM base as builder

# builder only dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive \
	apt-get install -y --no-install-recommends --no-install-suggests \
	liblzma-dev libsavvy-dev libshrinkwrap-dev \
	libsuperlu-dev libzstd-dev zlib1g-dev \
	r-cran-biocmanager r-bioc-biocversion r-cran-bh \
	r-cran-dplyr r-cran-qlcmatrix r-cran-r.utils \
	r-cran-rcpparmadillo r-cran-rcppeigen r-cran-rspectra \
	r-cran-rsqlite \
	&& rm -rf /var/lib/apt/lists/* \
	&& rm -rf /tmp/*

# pick up remaining bioconductor dependencies not packaged w/ ubuntu 
RUN Rscript -e 'BiocManager::install("SPAtest")' && \
	Rscript -e 'BiocManager::install("SKAT")' && \
	Rscript -e 'BiocManager::install("MetaSKAT")'

# copy over saige and remove thirdParty and configure script
# dependencies are already installed, so no need to run configure
COPY . SAIGE
RUN rm -rf SAIGE/thirdParty && \
	rm SAIGE/configure

RUN R CMD INSTALL SAIGE
# ------------------------------------------------------------------------------
# RUNTIME LAYER
FROM base
COPY --chmod=0555 --from=builder /SAIGE/extdata/*.R /usr/local/bin/
COPY --from=builder /usr/local/lib/R /usr/local/lib/R
# ------------------------------------------------------------------------------
