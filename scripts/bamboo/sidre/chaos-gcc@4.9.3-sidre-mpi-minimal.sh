#!/bin/bash
# This script builds a minimal set of the toolkit components designed for integration testing of sidre in host application codes.
# This build configuration is intended to mimick how a host code may build the sidre component.

# This version is parallel (mpi) and includes lumberjack and spio.
# It needs to be executed via an srun for the parallel tests to work.

if [[ $HOSTNAME == rz* ]]; then
    HC="host-configs/rzmerl-chaos_5_x86_64_ib-gcc@4.9.3.cmake"
else
    HC="host-configs/surface-chaos_5_x86_64_ib-gcc@4.9.3.cmake"
fi

BT="RelWithDebInfo"
BP="build-chaos-gcc@4.9.3-sidre-mpi-minimal"
IP="install-chaos-gcc@4.9.3-sidre-mpi-minimal"
COMP_OPT="-DENABLE_QUEST=OFF -DENABLE_SLAM=OFF -DENABLE_SHROUD=OFF"
BUILD_OPT="-DENABLE_DOCS=OFF -DENABLE_EXAMPLES=OFF -DENABLE_PYTHON=OFF -DENABLE_MPI=ON"


echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -ecc -hc $HC -bt $BT -bp $BP -ip $IP $COMP_OPT $BUILD_OPT    
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd $BP
	echo "Building..."
	echo "-----------------------------------------------------------------------"
	make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
	echo "-----------------------------------------------------------------------"

	echo "Run tests."
	echo "-----------------------------------------------------------------------"
	make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
	echo "-----------------------------------------------------------------------"
cd ..
