#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c clang@3.5.0 --buildtype RelWithDebInfo
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-clang@3.5.0-relwithdebinfo
    echo "Generating C/Fortran binding..."
    make generate
    if [ $? -ne 0 ]; then
        echo "Error: 'make generate' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Run tests"
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
    fi
    echo "-----------------------------------------------------------------------"
cd ..
