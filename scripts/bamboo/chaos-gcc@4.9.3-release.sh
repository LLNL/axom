#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Release --buildpath build-chaos-release --installpath install-chaos-release
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-release
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

    echo "Making docs..."
    echo "-----------------------------------------------------------------------"
    make VERBOSE=1 docs
    if [ $? -ne 0 ]; then
        echo "Error: 'make docs' failed"
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
cd ..

echo "Installing docs to web space..."
echo "-----------------------------------------------------------------------"
rm -rf /usr/global/web-pages/lc/www/toolkit/docs_old
if [ $? -ne 0 ]; then
    echo "Error: 'rm' failed"
    exit 1
fi

mv /usr/global/web-pages/lc/www/toolkit/docs /usr/global/web-pages/lc/www/toolkit/docs_old
if [ $? -ne 0 ]; then
    echo "Error: 'mv' failed"
    exit 1
fi

cp -R ./install-chaos-release/docs /usr/global/web-pages/lc/www/toolkit/
if [ $? -ne 0 ]; then
    echo "Error: 'cp' failed"
    exit 1
fi

chgrp -R toolkit /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chgrp' failed"
    exit 1
fi

chmod -R g+r+w+X /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
    exit 1
fi

chmod -R o+r+X /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
    exit 1
fi
echo "-----------------------------------------------------------------------"
