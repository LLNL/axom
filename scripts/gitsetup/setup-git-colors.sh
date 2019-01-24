#!/usr/bin/env bash

##
## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

##======================================================
#
# This script is used to set colors in git
#
##======================================================

# Project configuration instructions: NONE

read -ep 'Enable syntax highlighting of output from git commands? [Y/n] ' color &&
if test "$color" = "y" -o "$color" = "Y"; then
    git config --global color.ui true
else
    git config --global color.ui false
fi

