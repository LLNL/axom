#!/bin/bash
# 09-12-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug"
# 09-16-2016 chang28, build-and-test.sh "clang@3.5.0" "Debug" ""
# 09-19-2016 chang28, the decider has decided to have a configuration file call a main_script file, this is the main_script file, all environment variables are set up in the configuration file. 

#COMPILER="clang@3.5.0"
#COMPILER=$1

#if [[ $HOSTNAME == rz* ]]; then
#    HC="host-configs/rzmerl-chaos_5_x86_64_ib-${COMPILER}.cmake"
#else
#    HC="host-configs/surface-chaos_5_x86_64_ib-${COMPILER}.cmake"
#fi
#BT="Debug"
BT=$2
BP="atk_build"
IP="atk_install"
COMP_OPT=""
BUILD_OPT=$3
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

if ["$BUILD" = true]; then
    cd $BP
    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

if ["$TEST" = true]; then
    echo "Running tests..."
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test -j16"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

if ["$INSTALL_FILES" = true]; then
    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

f ["$INSTALL_FILES" = true]; then
    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"
fi

