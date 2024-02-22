#!/bin/bash

# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#BSUB -W 240
#BSUB -G wbronze
#BSUB -o b.out.cz.blueos.src.%J.txt
#BSUB -q pbatch
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_cz_blueos_src.sh
#

date
cd ..
./build_src.py
date
