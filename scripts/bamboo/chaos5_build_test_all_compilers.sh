#!/bin/bash
# 09-12-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug"
# 09-16-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug" ""
# 09-19-2016 chang28, the decider has decided to have a configuration file call a main_script file, this is the configuration file, all environment variables are set up here. chaos5_build_test_all_compilers.sh "Debug" ""


#BT="Debug"
BT=$1
BP="atk_build"
IP="atk_install"
COMP_OPT=""
BUILD_OPT=$2
OPTIONS="-ecc -hc $HC -bt $BT -bp $BP -ip $IP $COMP_OPT $BUILD_OPT"

BUILD=true
TEST=true
DOC=false
INSTALL_FILES=true
INSTALL_DOCS=false

TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/toolkit"
DOCS_DIR_OLD="${TOOLKIT_WEB_ROOT}/docs_old"
DOCS_DIR="${TOOLKIT_WEB_ROOT}/docs"

COMPILER="clang@3.5.0"
HC="host-configs/surface-chaos_5_x86_64_ib-${COMPILER}.cmake"
echo "Running "$COMPILER"
./main_script.sh
if [ $? -ne 0 ]; then
    echo "Error: 'calling ' $COMPILER '  failed"
    exit 1
fi
