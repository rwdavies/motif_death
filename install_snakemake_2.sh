## run interactively

cd /well/davies/users/dcc832/primates/bin/
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

##
## manually install conda to /well/davies/users/dcc832/primates/bin/miniconda3
## 

cd miniconda3/
. ./bin/activate
conda install -n base -c conda-forge mamba
conda activate base
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
snakemake --help


## activate and go!

conda install -n base -c conda-forge mamba
