# Automatic installation for HATCHet's preprocessing
: ex: set ft=markdown ;:<<'```shell' #

The following script allows the user to setup all the required requisites for running the preprocessing steps of HATCHet.

## Preliminary setup

First we define the `SETUP_HOME` directory as the current dictory, however any path can be alternatively chosen. Please update this path accoring to your requirements.

```shell
export SETUP_HOME=$(pwd)
:<<'```shell' # Ignore this line
```

We also ask the setup to terminate in case of errors and to print a trace of the execution by the following commands
```shell
set -e
set -o xtrace
PS4='[\t]'
:<<'```shell' # Ignore this line
```

## Downloading and installing SAMtools and BCFtools v1.6

We download and install SAMtools 1.6 within `${SETUP_HOME}/samtools-1.6`

```shell
wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 -P ${SETUP_HOME}/
tar -xvjf ${SETUP_HOME}/samtools-1.6.tar.bz2 -C ${SETUP_HOME}/ && rm -f ${SETUP_HOME}/samtools-1.6.tar.bz2
(cd ${SETUP_HOME}/samtools-1.6 && ./configure --prefix=${SETUP_HOME}/samtools-1.6/ && make -j && make -j install)
echo 'export PATH='${SETUP_HOME}'/samtools-1.6/bin/:${PATH}' > ${SETUP_HOME}/setup_hatchet.sh
:<<'```shell' # Ignore this line
```

We download and install BCFtools 1.6 within `${SETUP_HOME}/bcftools-1.6`

```shell
wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2 -P ${SETUP_HOME}/
tar -xvjf ${SETUP_HOME}/bcftools-1.6.tar.bz2 -C ${SETUP_HOME}/ && rm -f ${SETUP_HOME}/bcftools-1.6.tar.bz2
(cd ${SETUP_HOME}/bcftools-1.6 && ./configure --prefix=${SETUP_HOME}/bcftools-1.6/ && make -j && make -j install)
echo 'export PATH='${SETUP_HOME}'/bcftools-1.6/bin/:${PATH}'  >> ${SETUP_HOME}/setup_hatchet.sh
:<<'```shell' # Ignore this line
```

## Setting up the python environment


```shell
conda create -y -n hatchet python=3.8
# install HATCHet
conda activate hatchet
pip install -U pip
pip install -U setuptools
pip install .
# set up bioconda and install mosdepth
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda install mosdepth
echo "source ${SETUP_HOME}/miniconda/bin/activate hatchet" >> ${SETUP_HOME}/setup_hatchet.sh
:<<'```shell' # Ignore this line
```

## Using this setup for HATCHet

All the previous steps are one-time only, however the environments that we prepared should be activated in every new session before using HATCHet. This can be achieved by simply running the following command once per new session before running HATCHet, where `${SETUP_HOME}` is the path that we chose at the beginning:

```shell
source ${SETUP_HOME}/setup_hatchet.sh
exit $?
```
