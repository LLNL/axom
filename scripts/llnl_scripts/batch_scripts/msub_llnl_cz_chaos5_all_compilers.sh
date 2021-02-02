#!/bin/bash

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -l walltime=8:00:00
#MSUB -A wbronze
#MSUB -j oe
#MSUB -o m.out.cz.uberenv.chaos5.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_cz_chaos5_all_compilers.sh

date
cd ..
./build_tpls.py
date

