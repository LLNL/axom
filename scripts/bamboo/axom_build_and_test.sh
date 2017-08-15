#!/bin/bash
# 
# axom_build_and_test.sh "Debug" "" "clang@3.5.0" "develop" "cab" "chaos5" "/usr/workspace/wsrzd/axomdev/bamboo/rzmanta/build/clangcoral"
#
BUILD_TYPE=$1
BUILD_OPT=$2
COMPILER=$3
BRANCH=$4
HOST=$5
SYSTEM_TYPE=$6
WORKSPACE=$7

ssh $HOST "
   cd $WORSPACE && \
   echo 'Cloning...' && \
   echo '-----------------------------------------------------------------------' && \
   axom_codeCheckout.sh $BRANCH && \
   if [[ $SYSTEM_TYPE == 'chaos5' ] ||[ $SYSTEM_TYPE == 'toss3' ]]; then && \
      srun -N1 --exclusive -ppdebug ./scripts/bamboo/${SYSTEM_TYPE}_build_test_uno_compiler.sh  '$BUILD_TYPE' '$BUILD_OPT' 'COMPILER'; && \
   else && \
      ./scripts/bamboo/${SYSTEM_TYPE}_build_test_uno_compiler.sh '$BUILD_TYPE' '$BUILD_OPT'  '$COMPILER' && \
   fi && \
   "
if [ $? -ne 0 ]; then
   echo "Error: $SYSTEM_TYPE $HOST failed!"
   exit 1
fi

