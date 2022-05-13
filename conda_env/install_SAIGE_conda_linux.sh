#!/bin/bash
#

### RUN THIS SCRIPT IN A DIRECTORY IN WHICH YOU'RE HAPPY TO INSTALL SAIGE

# activate SAIGE env (use "source activate" if "conda activate" fails)
conda activate RSAIGE 2>/dev/null
if [ $? -ne 0 ]
then
   ACTPATH=`which activate 2>/dev/null` || ACTPATH=`whereis activate 2>/dev/null | awk '{print $2}'`
   DEACTPATH=`which deactivate 2>/dev/null` || DEACTPATH=`whereis deactivate 2>/dev/null | awk '{print $2}'`
   source $ACTPATH RSAIGE || { echo -e "Commands \"conda activate RSAIGE\" and \"source activate RSAIGE\" failed."; echo -e "Please check whether a) conda is installed and b) you have prepared the RSAIGE conda environment."; exit 1; }
fi

# set path for lib and compiler flags
FLAGPATH=`which python | sed 's|/bin/python$||'`

# set lib and compiler flags using python path
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"

# clone the SAIGE repo
src_branch=main
repo_src_url=https://github.com/saigegit/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url

# install MetaSKAT R package
echo "devtools::install_github(\"leeshawn/MetaSKAT\")" | R --slave

# install in R
R CMD INSTALL --library=$PWD/SAIGE SAIGE

# change the permissions of the executable files
chmod +x $PWD/SAIGE/extdata/step1_fitNULLGLMM.R
chmod +x $PWD/SAIGE/extdata/step2_SPAtests.R

# set PATH so that SAIGE scripts can be located on loading the conda env
# and set R library location so that SAIGE can be 
# then return PATH to normal on deactivating
mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
[ ! -s $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh ] && echo '#!/bin/sh' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
[ ! -s $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh ] && echo '#!/bin/sh' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo -e "export PATH_OLD=\$PATH\nexport PATH=\$PATH:$PWD/SAIGE/extdata\nexport R_LIBS_USER=$PWD/SAIGE" >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo -e "export PATH=\$PATH_OLD\nunset PATH_OLD\nunset R_LIBS_USER" >> $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

# deactivate the env
[ -z $DEACTPATH ] && conda deactivate || source $DEACTPATH
