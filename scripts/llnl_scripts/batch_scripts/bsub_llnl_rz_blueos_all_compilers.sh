#!/bin/bash

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#BSUB -W 240
#BSUB -o b.out.rz.blueos.all.compilers.%J.txt
#BSUB -q pdebug
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_all_compilers.sh
#

date
cd ..
./build_tpls.py
date

