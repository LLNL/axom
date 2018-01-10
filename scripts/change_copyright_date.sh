#!/usr/bin/env bash

##
## Copyright (c) 2017-2018-2018-2018, Lawrence Livermore National Security, LLC.
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

#=============================================================================
# Change the copyright date in all files that have the date.
# Just edit the 'grep' and 'sed' commands below to set what to search for
# and how to replace it.
#=============================================================================
#
# These are the commands you can use to replace the copyright date
# in all files.
#
# You may want to run each of these commands from the command line to
# make sure things are doing what you think they should be doing.
# This is why they are seperated into the steps; this could be more efficient.
# 

#=============================================================================
# First find all the files with old copyright dates
#=============================================================================
find . -type f ! -name \*.git\* -exec grep -l "This file is part of Axom" {} \; > files2change

#=============================================================================
# Replace the old copyright dates with new dates
#=============================================================================
for i in `cat files2change`
do
    echo $i
    cp $i $i.sed.bak
    sed "s/Copyright (c) 2017-2018/Copyright (c) 2017-2018/" $i.sed.bak > $i
done

#=============================================================================
# Remove the temporary files
#=============================================================================
find . -name \*.sed.bak -exec rm {} \;
rm files2change
