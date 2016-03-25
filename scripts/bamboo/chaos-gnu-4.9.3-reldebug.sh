#!/bin/sh
./scripts/config-build.py -c gcc@4.9.3 -bt RelWithDebInfo
cd build-chaos-gcc@4.9.3-relwithdebinfo
make -j16
make test ARGS="-T Test"

/usr/local/bin/ctest -T Test -R sidre

make install
make docs
