#!/bin/bash

echo "Configuring..."
echo "-----------------------------------------------------------------------"
./scripts/config-build.py -c gcc@4.9.3 --buildtype Release --buildpath build-chaos-release --installpath install-chaos-release
if [ $? -ne 0 ]; then
    echo "Error: config-build.py failed"
fi
echo "-----------------------------------------------------------------------"

cd build-chaos-release
    echo "Building..."
    echo "-----------------------------------------------------------------------"
    make -j16
    if [ $? -ne 0 ]; then
        echo "Error: 'make' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Running tests..."
    echo "-----------------------------------------------------------------------"
    make test ARGS="-T Test"
    if [ $? -ne 0 ]; then
        echo "Error: 'make test' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Making docs..."
    echo "-----------------------------------------------------------------------"
    make docs
    if [ $? -ne 0 ]; then
        echo "Error: 'make docs' failed"
    fi
    echo "-----------------------------------------------------------------------"

    echo "Installing files..."
    echo "-----------------------------------------------------------------------"
    make install
    if [ $? -ne 0 ]; then
        echo "Error: 'make install' failed"
    fi
    echo "-----------------------------------------------------------------------"
cd ..

echo "Installing docs to web space..."
echo "-----------------------------------------------------------------------"
rm -rf /usr/global/web-pages/lc/www/toolkit/docs_old
if [ $? -ne 0 ]; then
    echo "Error: 'rm' failed"
fi

mv /usr/global/web-pages/lc/www/toolkit/docs /usr/global/web-pages/lc/www/toolkit/docs_old
if [ $? -ne 0 ]; then
    echo "Error: 'mv' failed"
fi

cp -R ./install-chaos-release/docs /usr/global/web-pages/lc/www/toolkit/
if [ $? -ne 0 ]; then
    echo "Error: 'cp' failed"
fi

chgrp -R toolkit /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chgrp' failed"
fi

chmod -R g+r+w+X /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
fi

chmod -R o+r+X /usr/global/web-pages/lc/www/toolkit/docs/
if [ $? -ne 0 ]; then
    echo "Error: 'chmod' failed"
fi
echo "-----------------------------------------------------------------------"
