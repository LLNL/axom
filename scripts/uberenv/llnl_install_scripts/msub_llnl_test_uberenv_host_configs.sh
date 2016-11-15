#!/bin/bash
#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -l walltime=4:00:00
#MSUB -A wbronze
#MSUB -j oe
#MSUB -v UBERENV_PREFIX
#MSUB -o m.out.r.uberenv.test.atk.build.and.install.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
# env UBERENV_PREFIX={test path}  msub msub_llnl_test_uberenv_host_configs.sh

date
/usr/local/bin/python  llnl_lc_test_uberenv_host_configs.py
date

