#!/bin/bash
# 
# axom_launch.sh "Debug" "" "clang@3.5.0" 
#
echo  axom_launch.sh 1.0
set -ev
BUILD_TYPE=$1
BUILD_OPT=$2
COMPILER=$3

if [ $SYS_TYPE == chaos5* ]; then
   echo srun -N1 --exclusive -ppdebug ./scripts/bamboo/chaos5_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
elif [ $SYS_TYPE == toss3* ]; then
   echo srun -N1 --exclusive -ppdebug ./scripts/bamboo/toss3_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
else
   echo ./scripts/bamboo/bgo_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
fi

