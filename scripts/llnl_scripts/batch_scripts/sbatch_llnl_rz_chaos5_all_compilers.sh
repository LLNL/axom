#!/bin/bash

##
## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#SBATCH -N 1
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
cd ..
./build_tpls.py
date

