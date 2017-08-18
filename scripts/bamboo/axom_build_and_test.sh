#!/bin/bash
# 
# axom_build_and_test.sh "Debug" "" "clang@3.5.0" "cab" "/usr/workspace/wsb/axomdev/bamboo/cab/build/clang350/axom" 
#
echo  axom_build_and_test.sh 1.0
set -ev
BUILD_TYPE=$1
BUILD_OPT=$2
COMPILER=$3
HOST=$4
SYSTEM_TYPE=$5
WORKSPACE=$6


ssh $HOST "hostname; cd $WORKSPACE && \
   pwd && \
   ./scripts/bamboo/axom_launch.sh 'BUILD_TYPE' 'BUILD_OPT' 'COMPILER'"
