#!/bin/bash

# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#MSUB -l nodes=1
#MSUB -q pdebug
#MSUB -l walltime=4:00:00
#MSUB -A science
#MSUB -j oe
#MSUB -o m.out.rz.uberenv.bgqos.src.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_rz_bgq_src.sh
#

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use git-2.0.0
use cmake-3.8.2

cd ..
./build_src.py
date

