#!/bin/bash

COMPILER="gcc@4.9.3"

if [[ $HOSTNAME == rz* ]]; then
    HC="host-configs/rzmerl-chaos_5_x86_64_ib-${COMPILER}.cmake"
else
    HC="host-configs/surface-chaos_5_x86_64_ib-${COMPILER}.cmake"
fi

BT="Debug"
BP="build-chaos-${COMPILER}-${BT,,}"
IP="install-chaos-${COMPILER}-${BT,,}"
COMP_OPT=""
BUILD_OPT="-DENABLE_COVERAGE=TRUE"
OPTIONS="-ecc -hc $HC -bt $BT -bp $BP -ip $IP $COMP_OPT $BUILD_OPT"

echo "Configuring..."
echo "-----------------------------------------------------------------------"
echo "Options: $OPTIONS"
./scripts/config-build.py $OPTIONS
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd $BP
    echo "Generating C/Fortran binding..."
    make VERBOSE=1 generate
    if [ $? -ne 0 ]; then
        echo "Error: 'make generate' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Run tests and collect code coverage output."
    echo "-----------------------------------------------------------------------"
    make coverage ARGS="-T Test -E mpi"
    if [ $? -ne 0 ]; then
        echo "Error: 'make coverage' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    TOOLKIT_WEB_ROOT="/usr/global/web-pages/lc/www/toolkit"
    COVERAGE_DIR="${TOOLKIT_WEB_ROOT}/coverage"
    echo "Install code coverage report to web space..."
    echo "-----------------------------------------------------------------------"
    if [ -d ${COVERAGE_DIR} ]; then
        rm -rf ${COVERAGE_DIR}
        if [ $? -ne 0 ]; then
            echo "Error: 'rm' of coverage directory failed"
            #exit 1
        fi
    fi

    cp -R ./coverage ${TOOLKIT_WEB_ROOT}
    if [ $? -ne 0 ]; then
        echo "Error: 'cp' of coverage directory failed"
        exit 1
    fi
    
    chgrp -R toolkitd ${COVERAGE_DIR}
    if [ $? -ne 0 ]; then
        echo "Error: 'chgrp' on coverage directory failed"
        exit 1
    fi

    chmod -R ug+rwX ${COVERAGE_DIR}
    if [ $? -ne 0 ]; then
        echo "Error: 'chmod' on coverage directory failed"
        exit 1
    fi
    
    echo "-----------------------------------------------------------------------"

    echo "Running valgrind..."
    echo "-----------------------------------------------------------------------"
    make test ARGS=" -E mpi -D ExperimentalMemCheck --output-on-failure"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' with valgrind failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
cd ..
# TODO - Should we copy valgrind output to web area so team can access???
