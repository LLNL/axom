#!/bin/bash

##
## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
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
#  bsub  < bsub_llnl_cz_blueos_src.sh
#

date
cd ..
./build_src.py
date

