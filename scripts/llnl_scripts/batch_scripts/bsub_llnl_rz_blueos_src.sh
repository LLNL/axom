#!/bin/bash

# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#BSUB -n 8
#BSUB -W 240
#BSUB -G guests
#BSUB -x 
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_src.sh
#

date
cd ..
./build_src.py
date

