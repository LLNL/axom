#!/bin/bash

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#BSUB -N 1
#BSUB -W 240
#BSUB -G guests
#BSUB -o b.out.rz.blueos.all.compilers.%j.%N.txt
#BSUB -x 
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_all_compilers.sh
#

date
cd ..
./build_tpls.py
date

