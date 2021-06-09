#!/bin/bash

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#MSUB -l nodes=1:ppn=36
#MSUB -q pdebug
#MSUB -l walltime=4:00:00
#MSUB -j oe
#MSUB -o m.out.rz.uberenv.toss3.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_rz_toss3_all_compilers.sh

date
cd ..
./build_tpls.py
date

