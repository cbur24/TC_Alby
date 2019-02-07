#!/bin/bash
#PBS -P w85
#PBS -l walltime=2:00:00
#PBS -l mem=256GB
#PBS -l ncpus=32
#PBS -q normal
#PBS -m abe

#folder where scripts are stored

cd /g/data1a/w85/cb3058/wrf/TC_Alby

#load up an instance of python3 that has most of the libraries we need

module use /g/data/v10/public/modules/modulefiles/
module load dea-prod

# change the python-path to the location where you installed the packages that
# arent included in the standard instance of python.

export PYTHONPATH=/g/data1a/w85/cb3058/wrf/python_packages

# To install a python library in different location:
# 	pip install --target=/g/data1a/w85/cb3058/wrf/python_packages <packageName>


#Run the scripts. Make sure you've adjusted the user inputs sections.

python3 plot_wrfout_static.py >& plotting.log