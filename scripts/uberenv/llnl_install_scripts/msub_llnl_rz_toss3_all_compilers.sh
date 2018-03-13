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

#MSUB -l nodes=1:ppn=36
#MSUB -q pdebug
#MSUB -l walltime=8:00:00
#MSUB -j oe
#MSUB -o m.out.rz.uberenv.toss3.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_rz_toss3_all_compilers.sh

date
/usr/bin/python llnl_rz_uberenv_install_toss_3_x86_64_ib_all_compilers.py
date

