#!/bin/bash


COMPILER="intel@15.0.187"

if [[ $HOSTNAME == rz* ]]; then
    HC="host-configs/rzmerl-chaos_5_x86_64_ib-${COMPILER}.cmake"
else
    HC="host-configs/surface-chaos_5_x86_64_ib-${COMPILER}.cmake"
fi

BT="RelWithDebInfo"
BP="build-chaos-${COMPILER}-${BT,,}"
IP="install-chaos-${COMPILER}-${BT,,}"
COMP_OPT=""
BUILD_OPT="-DENABLE_CXX11=FALSE"



echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -ecc -hc $HC -bt $BT -bp $BP -ip $IP $COMP_OPT $BUILD_OPT    
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

    echo "Run tests"
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
cd..
