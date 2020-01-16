#!/bin/bash

# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

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

