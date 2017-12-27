#!/bin/bash

##
## Copyright (c) 2017, Lawrence Livermore National Security, LLC.
##
## Produced at the Lawrence Livermore National Laboratory.
##
## LLNL-CODE-741217
##
## All rights reserved.
##
## This file is part of Axom.
##
## For details about use and distribution, please read axom/LICENSE.
##

#BSUB -n 8
#BSUB -W 240
#BSUB -G guests
#BSUB -x 
#
# usage: 
#  cd {to directory with this script}
#  bsub  < bsub_llnl_rz_blueos_all_compilers.sh
#
# Note: Updates PYTHONPATH to include a local install of pyyaml
#       which is missing from LC's default python on blue_os

date
export PYTHONPATH=${PYTHONPATH}:/usr/gapps/axom/lib/python2.7/site-packages
cd ..
build_tpls.py
date

