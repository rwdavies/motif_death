#!/usr/bin/env bash

. activate

set -e

if [ ${BIN_DIR} == "" ]
then
    echo BIN_DIR not set
    exit 1
fi

mkdir -p ${BIN_DIR}
cd ${BIN_DIR} ## in global path

if [ ! -e virtualenv-15.1.0/virtualenv.py ]
then
    wget https://github.com/pypa/virtualenv/archive/15.1.0.tar.gz
    tar -xzvf 15.1.0.tar.gz
fi

export PATH=`pwd`/virtualenv-15.1.0/:${PATH}

rm -r -f snakemake
git clone https://github.com/snakemake/snakemake
cd snakemake
git checkout 55714744466391208f1355ab005990593f666501
virtualenv.py -p python3 .venv
source .venv/bin/activate
python setup.py install
