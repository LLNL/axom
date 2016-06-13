#/bin/bash
#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -l walltime=8:00:00
#MSUB -A wbronze
#MSUB -j oe
#MSUB -o m.out.r.uberenv.chaos5.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_rz_chaos5_all_compilers.sh

date
/usr/local/bin/python llnl_rz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py
date

