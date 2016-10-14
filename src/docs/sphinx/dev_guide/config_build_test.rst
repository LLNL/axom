.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

.. _configbuild-label:

======================================================
Configuration, Building, and Testing
======================================================

This section describes how to configure and build our code and how
to run tests and other code development tasks. See :ref:`repoclone-label`
for information about accessing the code.

.. note :: These sections need work...


.. _build-label:

--------------------------
Configuration and building
--------------------------

We use a CMake-based system to configure and build our code, called *BLT*
(see :ref:`tooleco-label` for more information). 

Python helper script
^^^^^^^^^^^^^^^^^^^^^

The easiest way to configure the code for compilation is to use the 
'config-build.py' python script in the 'scripts' directory; 
e.g.,::

   $ ./scripts/config-build.py -hc ./host-configd/surface-chaos_5_x86_64_ib-gcc@4.9.3.cmake

This script requires that you pass it a *host-config* file. The script runs 
CMake and passes it the host-config to initialize the CMake cache with the
configuration informarion contained in the file. See :ref:`hostconfig-label` 
for more information.

Running the script, as in the example above, creates two directories to hold
the build and install contents for the platform and compiler specified by the
host-config file - in this case, a CHAOS 5 platform with the GNU gcc 4.9.3
compiler. The The name 'surface' in the file name indicates the particular 
machine on which the host-config file was generated. Livermore Computing 
platforms are generally configured similarly so that the configuration will 
usually also work on other CHAOS 5 Linux platforms. 

To build the code and intall the header files, libraries, and documentation 
in the install directory, go into the build directory and run make; e.g.,::

   $ cd build-rzmerl-chaos_5_x86_64_ib-gcc@4.9.3-debug
   $ make
   $ make install

The python helper script accepts other arguments that allow you to specify
explicitly the build and install paths and build type. Following CMake 
conventions, we support three build types: 'Release', 'RelWithDebInfo', and 
'Debug'. To see the script options, run the script without any arguments; 
i.e.,::

   $ ./scripts/config-build.py 


Run CMake directly
^^^^^^^^^^^^^^^^^^^

You can also configure the code by running CMake directly and passing it 
the appropriate arguments. For example, to configure a release build with
the gcc compiler, you could do the following::

   $ mkdir build-gnu-release
   $ cd build-gnu-release
   $ cmake -DCMAKE_C_COMPILER=/usr/apps/gnu/4.9.3/bin/gcc \
     -DCMAKE_CXX_COMPILER=/usr/apps/gnu/4.9.3/bin/g++ \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gnu-release \
     ../src/
   $ make
   $ make install

.. note :: Actually, this will not work because the locations of several 
           third-party libraries must be provided. We should maintain a 
           list of dependencies required for the Toolkit components and
           make it clear what users need to provide to use what they need.


.. _hostconfig-label:

Host-config files
^^^^^^^^^^^^^^^^^^^

We use *host-config* files to track build configurations we support and 
maintain reproducibility. We maintain a collection of such files in the 
'host-configs' directory for platforms and compilers we support. 
When passed to CMake, using the '-C' option, a host-config file initializes 
the CMake cache with the configuration specified in the file. 

.. note :: Need to describe how the host config files get generated and how
           to generate new ones.



--------------------------
Make targets
--------------------------

Our system provides a variety of make targets to build individual Toolkit 
components, documentation, run tests, examples, etc. After running CMake 
(using either the python helper script or directly), you can see a listing of
all evailable targets by passing 'help' to make; i.e.,::

   $ make help

The name of each target should be sufficiently descriptive to indicate
what the target does. For example, to generate this developer guide, run the
following command::

   $ make dev_guide_docs


.. _testing-label:

--------------------------
Testing
--------------------------

Running tests
^^^^^^^^^^^^^^^

Describe how to run tests...

Adding tests
^^^^^^^^^^^^^^^

Describe how to add tests...



.. _tpl-label:

--------------------------
Third-party libraries
--------------------------

Describe how to run the scripts to install third-party libraries for 
testing different versions locally on a branch and for installing new
libraries for the team to use...

