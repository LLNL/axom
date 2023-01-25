#!/bin/bash
##############################################################################
# Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
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

echo HOST_CONFIG
echo $HOST_CONFIG

export BUILD_TYPE=${BUILD_TYPE:-Debug}


if [[ "$DO_BUILD" == "yes" ]] ; then
    echo "~~~~~~ RUNNING CMAKE ~~~~~~~~"
    or_die python3 ./config-build.py -hc /home/axom/axom/host-configs/docker/${HOST_CONFIG}.cmake -bt ${BUILD_TYPE} -DENABLE_GTEST_DEATH_TESTS=ON ${CMAKE_EXTRA_FLAGS}
    or_die cd build-$HOST_CONFIG-${BUILD_TYPE,,}
    echo "~~~~~~ BUILDING ~~~~~~~~"
    if [[ ${CMAKE_EXTRA_FLAGS} == *COVERAGE* ]] ; then
        or_die make -j 10
    else
        or_die make -j 10 VERBOSE=1
    fi
    if [[ "${DO_TEST}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING TESTS ~~~~~~~~"
        make CTEST_OUTPUT_ON_FAILURE=1 test ARGS='-T Test -VV -j8'
    fi
    if [[ "${DO_MEMCHECK}" == "yes" ]] ; then
        echo "~~~~~~ RUNNING MEMCHECK ~~~~~~~~"
        or_die ctest -T memcheck
    fi
fi

# Note: Azure pipelines requires read/write access for everyone between steps
find ./axom -type d -exec chmod 755 {} \;
find ./axom -type f -exec chmod 644 {} \;

if [[ "$DO_CLEAN" == "yes" ]] ; then
    echo "~~~~~~ CLEANING BUILD DIRECTORY ~~~~~~~~"
    make clean
fi

exit 0
