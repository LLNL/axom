#!/bin/bash
# 
# axom_build_and_test.sh "Debug" "" "clang@3.5.0" "cab" "chaos5" 
#
BUILD_TYPE=$1
BUILD_OPT=$2
COMPILER=$3
HOST=$4
SYSTEM_TYPE=$5

if [ $SYSTEM_TYPE == "chaos5" ] || [ $SYSTEM_TYPE == "toss3" ]; then
   SCRIPT="srun -N1 --exclusive -ppdebug ./scripts/bamboo/${SYSTEM_TYPE}_build_test_uno_compiler.sh";
else
   SCRIPT="./scripts/bamboo/${SYSTEM_TYPE}_build_test_uno_compiler.sh ";
fi
echo $SCRIPT

cur_dir=`pwd`
ssh $HOST "hostname; cd $cur_dir && \
   pwd && \
   $SCRIPT '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER'"

if [ $? -ne 0 ]; then
   echo "Error: $SYSTEM_TYPE $HOST failed!"
   exit 1
fi

