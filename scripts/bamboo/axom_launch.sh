#!/bin/bash
# 
# axom_launch.sh "Debug" "" "clang@3.5.0" 
#
echo  axom_launch.sh 1.0
set -ev
BUILD_TYPE=$1
BUILD_OPT=$2
COMPILER=$3

pwd
echo $SYS_TYPE
if [ $SYS_TYPE == "chaos_5_x86_64_ib" ]; then
   srun -N1 --exclusive -ppdebug ./scripts/bamboo/chaos5_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
elif [ $SYS_TYPE == "toss_3_x86_64_ib" ]; then
   srun -N1 --exclusive -ppdebug ./scripts/bamboo/toss3_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
elif [ $SYS_TYPE == "bgqos_0" ]; then
   ./scripts/bamboo/bgq_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
elif [ $SYS_TYPE == "blueos_3_ppc64le_ib" ]; then
    ./scripts/bamboo/blueos_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'
fi

