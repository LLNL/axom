#!/bin/bash
#SBATCH -N 8
#SBATCH -J axom_uberenv_rz_chaos5
#SBATCH -t 8:00:00
#SBATCH -p pdebug
#SBATCH -A wbronze
#SBATCH --exclusive
#SBATCH -o m.out.sbatch.rz.uberenv.chaos5.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  sbatch sbatch_llnl_rz_chaos5_all_compilers.sh

date
/usr/local/bin/python llnl_rz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py
date

