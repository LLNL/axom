#!/bin/bash

set -ev 

WORKSPACE=`pwd`/../..
cd $WORKSPACE

echo "~~~~~~~~~~~~ Building tpls ~~~~~~~~~~~"
cd ./scripts/uberenv
/usr/bin/python ./uberenv.py --spec="^conduit@master" --prefix="uberenv_libs" # uberenv_libs is the default dir name

echo "~~~~~~~~~~~~ Building axom ~~~~~~~~~~~"
cd ../..
./config-build.py -hc ./scripts/uberenv/uberenv_libs/*.llnl.gov-darwin-x86_64-clang@9.0.0.cmake -bp ./build
cd build
make

echo "~~~~~~~~~~~~ Testing axom ~~~~~~~~~~~"
make test

