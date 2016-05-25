#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -hc host-configs/surface-chaos_5_x86_64_ib-gcc\@4.9.3.cmake   --buildtype Debug --buildpath atk_build --installpath atk_install
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd atk_build
    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Running tests..."
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
        exit 1
    fi
    echo "-----------------------------------------------------------------------"

