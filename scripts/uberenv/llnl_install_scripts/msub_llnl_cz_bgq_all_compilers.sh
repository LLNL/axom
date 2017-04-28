#!/bin/bash
#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -l walltime=8:00:00
#MSUB -A ddcwork
#MSUB -j oe
#MSUB -o m.out.cz.uberenv.bgqos.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_cz_bgq_all_compilers.sh

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use cmake-3.1
python ./llnl_cz_uberenv_install_bgqos_0_all_compilers.py
date

