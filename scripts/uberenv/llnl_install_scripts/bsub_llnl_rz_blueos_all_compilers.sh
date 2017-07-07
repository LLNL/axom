#!/bin/bash
#BSUB -n 1
#BSUB -q pdebug
#BSUB -w 240
#BSUB -A guests
#BSUB -x 
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_all_compilers.sh

date
python ./llnl_rz_uberenv_install_blueos_3_ppc64le_ib_all_compilers.py
date

