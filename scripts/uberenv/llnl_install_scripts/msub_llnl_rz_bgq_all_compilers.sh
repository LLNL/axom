#!/bin/bash

##
## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## LLNL-CODE-741217
##
## All rights reserved.
##
## This file is part of Axom.
##
## For details about use and distribution, please read axom/LICENSE.
##

#MSUB -l nodes=1
#MSUB -q pdebug
#MSUB -l walltime=8:00:00
#MSUB -A science
#MSUB -j oe
#MSUB -o m.out.rz.uberenv.bgqos.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_rz_bgq_all_compilers.sh

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use git-2.0.0
use cmake-3.8.2
python ./llnl_rz_uberenv_install_bgqos_0_all_compilers.py
date

