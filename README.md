# SAIGE-GPU
This GitHub repository branch contains the GPU version of SAIGE version 1.3.3.
We have modified the original version of SAIGE step 1, by building the full GRM through distribution to multiple GPUs.

This version was used to run the VA MVP gwPheWAS on the DOE ORNL OLCF Summit HPC.

## Install via Conda

### Prerequisites

**Conda**: Conda should be available on Polaris.

### Installation Steps

#### Step 1: Create directory and clone the SAIGE-DOE GitHub repository

Create and activate a new Conda environment for SAIGE:

```bash
mkdir SAIGE-GPU
cd SAIGE-GPU
git clone https://github.com/saigegit/SAIGE.git
cd SAIGE/
git pull origin SAIGE-GPU-1.3.3
```

#### Step 2: Create a Conda Environment

Create and activate a new Conda environment for SAIGE. We will be using the existing YML file from the repository.
In this example I am naming my environment `RSAIGE_GPU`, you can replace with your desired name. 
In addition, due to the size of the conda environment and packages, I will be installing the environment in a different mount where I have more space. You can create your environment either in your home directory if you have enough space, or your project space.

```bash
conda env create  --file=./conda_env/environment-RSAIGE.yml -p /home/user/conda_envs/RSAIGE_GPU
conda activate /home/user/conda_envs/RSAIGE_GPU
```

#### Step 3: Install Dependencies
SAIGE requires several packages and libraries. Install these with the following commands:


```bash
pip3 install cget click
conda install cuda -c nvidia/label/cuda-11.4.3
```

You need to install openMPI version `4.1.5`. This is difficult to perform within Conda, so we will install separately, but then include it in our Conda environment. To do so, we will first deactivate the `conda` environment and install `openmpi-4.1.5`:

```bash
conda deactivate
cd /home/user/conda_envs/pkgs
wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz
tar -xzf openmpi-4.1.5.tar.gz
cd openmpi-4.1.5
./configure --prefix=/home/user/conda_envs/RSAIGE_GPU/opt/openmpi
make -j4
make install
export PATH=/home/user/conda_envs/RSAIGE_GPU/opt/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/user/conda_envs/RSAIGE_GPU/opt/openmpi/lib

mkdir -p $CONDA_PREFIX/etc/conda/activate.d
mkdir -p $CONDA_PREFIX/etc/conda/deactivate.d
echo 'export PATH=$CONDA_PREFIX/opt/openmpi/bin:$PATH' > $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
echo 'export LD_LIBRARY_PATH=$CONDA_PREFIX/opt/openmpi/lib:$LD_LIBRARY_PATH' >> $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
```

Now you can activate the `conda` environment once again. 
Install other required libraries, such as pbdMPI, savvy, superlu:

```bash
conda activate /home/user/conda_envs/RSAIGE_GPU
Rscript -e 'install.packages("pbdMPI", repos=c("https://cloud.r-project.org"))'
conda install -c conda-forge -c bioconda savvy
conda install conda-forge::superlu
```

#### Step 4: Compile SAIGE
To compile SAIGE-GPU, first clean any previous builds and then run make:

```bash
cd ./SAIGE-GPU
R CMD INSTALL SAIGE-DOE
```

If you encounter linking errors, ensure that the PKG_LIBS line in the Makevars file correctly references the MPI library.

#### Step 5: Verify Installation
Check if the installation for step 1 was successful by running the following commands where the output should be the list of parameter options:

```bash
path_to_saige=./SAIGE-GPU/SAIGE-DOE
Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R --help
```

If the help information is displayed for each command, the installation is complete.

You can also run a test with the provided test input files. You can replace `mpirun -n 4` in the command below with the appropriate number of GPUs you have available.

```bash
# access a node with GPUs
# once the node is provided
conda activate /home/user/conda_envs/RSAIGE_GPU

path_to_saige=~/SAIGE-GPU/SAIGE-DOE
mpirun -n 4 Rscript $path_to_saige/extdata/step1_fitNULLGLMM.R \
--plinkFile=$path_to_saige/extdata/input/plinkforGRM_1000samples_10kMarkers \
--phenoFile=$path_to_saige/extdata/input/pheno_1000samples.txt \
--invNormalize=FALSE \
--phenoCol=y \
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=./GPU_step1_output \
--minMAFforGRM 0.01 \
--LOCO  \
 --IsOverwriteVarianceRatioFile=TRUE \
--nThreads=1
```

You should see a succesful run where all GPUs are used. The log should provide the IDs of the GPUs used.

Currently, Step 2 is not optimized to use GPUs and but can run multiple phenotypes in parallel using nThreads. 
You can also test step 2 for a single phenotype by running the following:

```
 Rscript step2_SPAtests.R        \
     --bgenFile=./input/genotype_100markers.bgen    \
     --bgenFileIndex=./input/genotype_100markers.bgen.bgi \
     --SAIGEOutputFile=./genotype_100markers_bgen_Test_out_GRMforStep1.txt \
     --chrom=1 \
     --AlleleOrder=ref-first \
     --minMAF=0 \
     --minMAC=0.5 \
     --sampleFile=./input/samplelist.txt \
     --GMMATmodelFile=./GPU_step1_output.rda \
     --varianceRatioFile=./GPU_step1_output.varianceRatio.txt   \
     --LOCO=FALSE       \
     --is_fastTest=TRUE \
     --nThreas 2

```

## Run on HPC Systems Using Singularity Container

It is not necessary to build SAIGE-GPU from source code. Instead you can pull the latest Docker container and run a small example from the SAIGE package itself:

```
singularity build saige-doe.sif docker://tnnandi/saige-doe:2

mpirun -n <number of gpus>
   singularityexec --bind /SAIGE_container/SAIGE-DOE/extdata saige_1.1.9.sif /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
   --plinkFile=$ARRAYS_DIR/100k_arrays \
   --phenoFile=$INPUT_DIR/phenotypes.tsv \
   --invNormalize=FALSE \
   --phenoCol=case \
   --covarColList=sex_at_birth_Male,age,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
   --qCovarColList=sex_at_birth_Male \
   --sampleIDColinphenoFile=IID \
   --traitType=binary \
   --outputPrefix=$TMPDIR/GPU_step1_output \
   --minMAFforGRM 0.01 \
   --LOCO FALSE \
   --IsOverwriteVarianceRatioFile=TRUE \
   --nThreads=1; gsutil -m cp $TMPDIR/GPU_step1_output* $OUT_DIR/
```



## Run on DNANexus for UK Biobank

This guide provides step-by-step instructions to execute **SAIGE-GPU Step 1** for genome-wide association studies (GWAS) using GPU-accelerated tools for the DNAnexus platform for UKBB data. It includes memory calculation, dataset preparation, and command-line execution details.

---

### 1. Calculate Memory Requirements

To estimate the memory required for Step 1:

- **Formula**:  
  Memory (GB) = `4 * M * N / 1e+9`  
  Where:
  - `M`: Number of genetic markers (variants).
  - `N`: Number of samples (participants).

- **Example Calculation**:  
  For a GWAS with 112,871 variants and 435,212 samples:  
  ```plaintext
  (4 * 435,212 * 112,871) / 1e+9 ≈ 196.49 GB
  ```
- **Note**: Ensure an additional 20 GB of memory for auxiliary requirements.
- **Check GPU Memory**:
  Use the nvidia-smi command to assess available GPU memory.
  
### 2. Launch a Cloud Workstation Job

Prepare the following inputs for Step 1 on DNAnexus or another cloud platform:

- **Required Files**:
  - SAIGE-GPU image or tarball.
  - PLINK files:
      - .bed
      - .bim
      - .fam
  - Phenotype-covariate table.
 
### 3. Install Docker

Install Docker on your system by following [Docker’s installation guide for Ubuntu](https://docs.docker.com/engine/install/ubuntu/).

### 4. OPTIONAL: Prune Datasets for Memory Requirements

If your dataset exceeds memory limits, prune it as follows:

- **Download PLINK**:

```bash
  wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20241022.zip
  unzip plink_linux_x86_64_20241022.zip
```

- **Randomly select variants**:

```bash
  shuf -n [number of variants] [your plink dataset prefix].bim > random_variants.txt
```

- **Generate New PLINK Files**:

```bash
  ./plink --bfile [your PLINK dataset prefix] --extract random_variants.txt --make-bed --out [new PLINK dataset prefix]
```

### 5. Prepare Files and Load Docker Image

- **Pull Docker Image**:

```bash
    docker pull /tnnandi/saige-doe:2
```

- **Save image for faster use**:

```bash
   docker save tnnandi/saige-doe:2 | gzip > saige_gpu_image.tar.gz
```

- **Load Docker Image**:

```bash
   docker load -i saige_gpu_image.tar.gz
```

- **Organize Data**:

```bash
  mkdir data experiment
  mv [your PLINK dataset prefix].* data/
  mv [your pheno-covar table] data/
```

### 6. Run Step 1

Replace [ ] with your specific dataset information and run the following command:

```bash
  docker run --rm \
           -v /home/dnanexus/data:/SAIGE_container/data \
           -v $PWD/experiment:/SAIGE_container/output \
           --gpus device=all \
           saige_gpu_image \
           mpirun -n [number of gpus] --allow-run-as-root \
           /opt/conda/lib/R/bin/Rscript /SAIGE_container/SAIGE-DOE/extdata/step1_fitNULLGLMM.R \
           --plinkFile=/SAIGE_container/data/[your PLINK dataset prefix] \
           --phenoFile=/SAIGE_container/data/[your pheno-covar table] \
           --phenoCol=[Phenotype column name] \
           --covarColList=[Covariates list] \
           --sampleIDColinphenoFile=[Sample ID column name] \
           --traitType=binary \
           --outputPrefix=/SAIGE_container/output/step1 \
           --nThreads=1 \
           --minMAFforGRM [Minimum MAF value] \
           --maxiterPCG [Max PCG iterations] \
           --maxiter [Max iterations] \
           --LOCO FALSE \
           --IsOverwriteVarianceRatioFile=TRUE > step1.log
```

### 7. Monitoring and Troubleshooting

- **Run in `tmux` for monitoring**:
  - Monitor GPU usage:
    ```bash
       nvidia-smi
    ```
    
  - Check SAIGE log output:
    ```bash
       tail -f step1.log
    ```




Please view the [SAIGE-GPU pages](https://exascale-genomics.github.io/SAIGE-GPU/index.html#/) for a more detailed description.


#######
HERE we are maintaining an newly improved stable version of SAIGE and SAIGE-GENE+. 
Please find the https://saigegit.github.io/SAIGE-doc/ for documentation.




