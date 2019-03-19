#!/usr/bin/env bash

###############################################################################
# Copyright (c) 2017-19, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#       
# All rights reserved.
#   
# This file is part of Axom.
#   
# For details about use and distribution, please read axom/LICENSE.
#   
###############################################################################

TAR_CMD=`which tar`
VERSION=`git describe --tags`

git archive --prefix=Axom-${VERSION}/ -o Axom-${VERSION}.tar HEAD 2> /dev/null

echo "Running git archive submodules..."

p=`pwd` && (echo .; git submodule foreach --recursive) | while read entering path; do
    temp="${path%\'}";
    temp="${temp#\'}";
    path=$temp;
    [ "$path" = "" ] && continue;
    (cd $path && git archive --prefix=Axom-${VERSION}/$path/ HEAD > $p/tmp.tar && ${TAR_CMD} --concatenate --file=$p/Axom-${VERSION}.tar $p/tmp.tar && rm $p/tmp.tar);
done

gzip Axom-${VERSION}.tar
