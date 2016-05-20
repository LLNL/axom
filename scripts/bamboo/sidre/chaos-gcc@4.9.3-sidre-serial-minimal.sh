#!/bin/bash
# This script builds a minimal set of the toolkit components designed for integration testing of sidre in host application codes.
# This build configuration is intended to mimick how a host code may build the sidre component.

# This version is serial (no mpi).

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype RelWithDebInfo -DBUILD_DOCS=OFF -DBUILD_EXAMPLES=OFF -DENABLE_QUEST=OFF -DENABLE_SLAM=OFF -DENABLE_SHROUD=OFF -DENABLE_PYTHON=OFF -DENABLE_MPI=OFF
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-gcc@4.9.3-relwithdebinfo
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
