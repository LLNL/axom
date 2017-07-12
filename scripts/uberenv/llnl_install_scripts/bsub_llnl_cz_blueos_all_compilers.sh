#!/bin/bash
#BSUB -n 1
#BSUB -W 240
#BSUB -G guests
#BSUB -x 
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_all_compilers.sh

date
python ./llnl_cz_uberenv_install_blueos_3_ppc64le_ib_all_compilers.py
date

