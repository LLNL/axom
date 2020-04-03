#!/bin/bash
##############################################################################
# Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
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

echo HOST_CONFIG
echo $HOST_CONFIG

if [[ "$DO_BUILD" == "yes" ]] ; then
    or_die ./config-build.py -hc /home/axom/axom/host-configs/docker/${HOST_CONFIG}.cmake -DENABLE_GTEST_DEATH_TESTS=ON
    or_die cd build-$HOST_CONFIG-debug
    if [[ ${CMAKE_EXTRA_FLAGS} == *COVERAGE* ]] ; then
        or_die make -j 8
    else
        or_die make -j 8 VERBOSE=1
    fi
    if [[ "${DO_TEST}" == "yes" ]] ; then
        make CTEST_OUTPUT_ON_FAILURE=1 test ARGS='-T Test -VV -j8'
    fi
    if [[ "${DO_MEMCHECK}" == "yes" ]] ; then
      or_die ctest -T memcheck
    fi
fi

find ./axom/sidre -type d -exec chmod 755 {} \;
find ./axom/sidre -type f -exec chmod 644 {} \;

exit 0
