#!/bin/bash
# 09-12-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug"
# 09-16-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug" ""
# 09-19-2016 chang28, the decider has decided to have a configuration file call a main_script file, this is the configuration file, all environment variables are set up here. chaos5_build_test_all_compilers.sh "Debug" ""
# 09-20-2016 chang28,  chaos5_build_test_uno_compiler.sh "Debug" "" "clang@3.5.0"

echo chaos5_build_and_test_uno_compiler version 0.9.9
#BT="Debug"
BUILD_TYPE=$1
BUILD_PATH="atk_build"
INSTALL_PATH="atk_install"
COMP_OPT=""
BUILD_OPT=$2

BUILD=true
TEST=true
DOC=false
INSTALL_FILES=true
INSTALL_DOCS=false

JOBS=16

TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/toolkit"
DOCS_DIR_OLD="${TOOLKIT_WEB_ROOT}/docs_old"
DOCS_DIR="${TOOLKIT_WEB_ROOT}/docs"

COMPILER=$3
if [[ $HOSTNAME == rz* ]]; then
    HOST_CONFIGURATION="host-configs/rzmerl-chaos_5_x86_64_ib-${COMPILER}.cmake"
else
    HOST_CONFIGURATION="host-configs/cab-chaos_5_x86_64_ib-${COMPILER}.cmake"
fi

OPTIONS="-ecc -hc $HOST_CONFIGURATION -bt $BUILD_TYPE -bp $BUILD_PATH -ip $INSTALL_PATH $COMP_OPT $BUILD_OPT"
echo Running $COMPILER
. ./scripts/bamboo/main_script.sh
if [ $? -ne 0 ]; then
    echo Error: calling  $COMPILER  failed
    exit 1
fi
