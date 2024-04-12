#!/bin/bash
##############################################################################
# Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level LICENSE file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)
##############################################################################

set -x

function or_die () {
    "$@"
    local status=$?
    if [[ $status != 0 ]] ; then
        echo ERROR $status command: $@
        exit $status
    fi
}

or_die cd axom
git submodule init
git submodule update

echo "~~~~ helpful info ~~~~"
echo "USER="`id -u -n`
echo "PWD="`pwd`
echo "~~~~~~~~~~~~~~~~~~~~~~"

echo "~~~~~~~~~ls -al /~~~~~~~~~~"
ls -al /
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home~~~~~~~~~~"
ls -al /home
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home/axom~~~~~~~~~~"
ls -al /home/axom
echo "~~~~~~~~~~~~~~~~~~~~~~"
echo "~~~~~~~~~ls -al /home/axom/axom~~~~~~~~~~"
ls -al /home/axom/axom
echo "~~~~~~~~~~~~~~~~~~~~~~"


echo "~~~~~~ RUNNING SPACK ~~~~~~~~"
UPSTREAM=/home/axom/axom_tpls/spack
ls -al $UPSTREAM
python3 /home/axom/axom/scripts/uberenv/uberenv.py \
    --upstream=$UPSTREAM \
    --spack-env-file=/home/axom/axom/scripts/spack/configs/linux_ubuntu_20/spack.yaml \
    --spack-build-mode=dev-build
    --package-final-phase=install \
    --run_tests

exit 0
