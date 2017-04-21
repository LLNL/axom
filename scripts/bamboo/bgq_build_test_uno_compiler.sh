#!/bin/bash
# 10-20-2016 chang28, 
# 10-20-2016 chang28, bgq_build_test_uno_compiler.sh "Debug" "" "gcc@4.7.2"
# 01-13-2017 chang28, turn off testing
# 03-22-2017 weiss27, Need to run with two configuration options on BG/Q and compile on login nodes
#                     bgq_build_test_uno_compiler.sh "Debug" "" gcc@4.7.2

echo bgq version 1.0.0
#BT="Debug"
BUILD_TYPE=$1
BUILD_PATH="atk_build"
INSTALL_PATH="atk_install"
COMP_OPT=""
BUILD_OPT=$2

BUILD=true
TEST=false
DOC=true
INSTALL_FILES=true
INSTALL_DOCS=false

# Explicitly set a CTest that works on the backend nodes 
CTEST_EXE="/usr/global/tools/CMake/bgqos_0/cmake-3.0-bgq-experimental/bin/ctest"  
RUN_BGQ_TEST=true

JOBS=16

TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/axom"
DOCS_DIR_OLD="${TOOLKIT_WEB_ROOT}/docs_old"
DOCS_DIR="${TOOLKIT_WEB_ROOT}/docs"

COMPILER=$3
if [[ $HOSTNAME == rz* ]]; then
    HOST_CONFIGURATION="host-configs/rzuseqlac-bgqos_0-${COMPILER}.cmake"
    PARTITION="pdebug"
else
    HOST_CONFIGURATION="host-configs/vulcanlac-bgqos_0-${COMPILER}.cmake"
    PARTITION="psmall"
fi

OPTIONS="-ecc -hc $HOST_CONFIGURATION -bt $BUILD_TYPE -bp $BUILD_PATH -ip $INSTALL_PATH $COMP_OPT $BUILD_OPT"
echo Running $COMPILER
. ./scripts/bamboo/main_script.sh
if [ $? -ne 0 ]; then
    echo Error: calling  $COMPILER  failed
    exit 1
fi

# Now run the tests using a ctest executable that we know to work on BG/Q
if [ "$RUN_BGQ_TEST" = true ]; then
    cd $BUILD_PATH
    pwd
    echo "Running tests env CTEST_OUTPUT_ON_FAILURE=1 ..."
    export CTEST_OUTPUT_ON_FAILURE=1
    echo "-----------------------------------------------------------------------"
    salloc -N4 -p$PARTITION $CTEST_EXE -T Test -j$JOBS
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    cd ..
    echo "-----------------------------------------------------------------------"
fi

