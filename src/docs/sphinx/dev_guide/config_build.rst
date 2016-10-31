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
Configure and Build
======================================================

This section describes more advanced CS Toolkit build and configuration 
tasks that all developers should be aware of and be able to perform. 
More basic build instructions can be found in the User Quick Start Guide
**add link to that guide**. See :ref:`repoclone-label` for information 
about accessing the code.

We use a CMake-based system to configure and build our code, called *BLT*
(see :ref:`tooleco-label` for more information). 

---------------------
Python Helper Script
---------------------

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
compiler. The name 'surface' in the file name indicates the particular 
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

-----------------------
Running CMake Directly
-----------------------

You can also configure the code by running CMake directly and passing it 
the appropriate arguments. For example, to configure, build and install 
a release build with the gcc compiler, you could pass a host-config file 
CMake::

   $ mkdir build-gnu-release
   $ cd build-gnu-release
   $ cmake -C ./host-configd/surface-chaos_5_x86_64_ib-gcc@4.9.3.cmake \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gnu-release \
     ../src/
   $ make
   $ make install

You can also run CMake by explicitly passing all options you need. Here is 
a summary of commonly used CMake options:

.. note:: **Fill this in...** 


.. _hostconfig-label:

------------------
Host-config Files
------------------

We use *host-config* files to track build configurations we support and 
maintain reproducibility. We maintain a collection of such files in the 
'host-configs' directory for platforms and compilers we support. 
When passed to CMake, using the '-C' option, a host-config file initializes 
the CMake cache with the configuration specified in the file. 

.. note :: Need to describe how the host config files get generated and how
           to generate new ones.


--------------------------
Make Targets
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

.. note :: Add a table that provides an overview of our make targets.


.. _tpl-label:

--------------------------
Third-party Libraries
--------------------------

Describe how to run the scripts to install third-party libraries for 
testing different versions locally on a branch and for installing new
libraries for the team to use...

Building and installing TPLs for all compilers on LC CHAOS platforms (CZ)::

   $ python ./scripts/uberenv/llnl_install_scripts/llnl_cz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py

Questions we need to answer include:

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs with for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    experimentation?
  * Others?

.. note :: Pull in content from ../web/build_system/thirdparty_deps.rst ...
           fill in gaps and make sure it it up-to-date...
           

