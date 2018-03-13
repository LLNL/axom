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
cd ..
./build_src.py
date

