#!/bin/bash

# Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#MSUB -l nodes=1
#MSUB -q pbatch
#MSUB -l walltime=4:00:00
#MSUB -A wbronze
#MSUB -j oe
#MSUB -v UBERENV_PREFIX
#MSUB -o m.out.r.uberenv.test.atk.build.and.install.%j.%N.txt
#
# usage: 
#  cd {to directory with this script}
# env UBERENV_PREFIX={test path}  msub msub_llnl_test_uberenv_host_configs.sh

date
cd ..
./build_src.py
date

