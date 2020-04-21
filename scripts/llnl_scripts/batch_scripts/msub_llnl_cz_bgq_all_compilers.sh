#!/bin/bash

# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#MSUB -l nodes=1
#MSUB -q psmall
#MSUB -l walltime=8:00:00
#MSUB -A ddcwork
#MSUB -j oe
#MSUB -o m.out.cz.uberenv.bgqos.all.compilers.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
#  msub -d `pwd` msub_llnl_cz_bgq_all_compilers.sh
#

date
. /usr/local/tools/dotkit/init.sh
use python-2.7.3
use git-2.0.0
use cmake-3.8.2

cd ..
./build_tpls.py
date

