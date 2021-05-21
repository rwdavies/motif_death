. activate ## includes newer python, virtualenv
mkdir ${ANALYSIS_DIR}bin
cd ${ANALYSIS_DIR}bin
set -e 
git clone https://bitbucket.org/snakemake/snakemake.git
cd snakemake
virtualenv -p python3 .venv
source .venv/bin/activate
python setup.py install
