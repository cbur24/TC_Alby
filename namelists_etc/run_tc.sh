#!/bin/bash
#PBS -P w85
#PBS -l walltime=1:00:00
#PBS -l mem=64GB
#PBS -l ncpus=16
#PBS -q express

cd /short/w85/cb3058/WRFV_3.7.1/WRFV3/test/em_real

module load dot
module load intel-fc/12.1.9.293
module load intel-cc/12.1.9.293
module load openmpi/1.6.3
module load netcdf/4.3.3.1
export JASPERINC=/usr/include
export JASPERLIB=/usr/lib64
export WRFIO_NCD_LARGE_FILE_SUPPORT=1

./tc.exe >& TCoutput.log

