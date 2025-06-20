FROM ubuntu:24.04

LABEL maintainer="wzhou@broadinstitute.org"

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV OMP_NUM_THREADS=1
# Force step_2 to use 1 single thread. More threads are ineffective
ENV PIXI_NO_PATH_UPDATE=1

RUN apt-get update -qq && \
    apt-get -y install curl time \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . .

RUN curl -fsSL https://pixi.sh/install.sh | sh && \
    mv /root/.pixi/bin/pixi /bin && pixi install && \
    rm -rf /root/.cache && \
    pixi run Rscript -e 'install.packages("lintools", repos="https://cloud.r-project.org")'

RUN curl -L https://github.com/chrchang/plink-ng/archive/refs/tags/v2.0.0-a.6.16.tar.gz | tar -zx && \
    mv plink-ng-2.0.0-a.6.16 plink-ng && \
    pixi run x86_64-conda-linux-gnu-cc -std=c++14 -fPIC -O3 -I.pixi/envs/default/include -L.pixi/envs/default/lib -o plink2_includes.a plink-ng/2.0/include/*.cc -shared -lz -lzstd -lpthread -lm -ldeflate && \
    mv plink2_includes.a .pixi/envs/default/lib

RUN pixi run R CMD INSTALL .

RUN mv extdata/step1_fitNULLGLMM.R extdata/step2_SPAtests.R extdata/step3_LDmat.R extdata/createSparseGRM.R /usr/local/bin/ && \
   chmod a+x /usr/local/bin/step1_fitNULLGLMM.R /usr/local/bin/step2_SPAtests.R /usr/local/bin/step3_LDmat.R /usr/local/bin/createSparseGRM.R && \
   createSparseGRM.R --help && \
   step1_fitNULLGLMM.R --help && \
   step2_SPAtests.R --help && \
   step3_LDmat.R --help