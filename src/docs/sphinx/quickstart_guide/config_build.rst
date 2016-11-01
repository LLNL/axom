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

======================================================
Configuration and Building
======================================================

This section provides information about configuring and building
the CS Toolkit software after you have obtained a copy of the repo.
The main steps for using the Toolkit are:

  #. Configure, build, and install third-party libraries (TPLs) that the Toolkit depends on.
  #. Build and install the CS Toolkit and its libraries.
  #. Build and link your application with the Toolkit installation.

Depending on how your team uses the Toolkit, some of these steps, such as
installing the Toolkit TPLs and the Toolkit itself, may need to be done 
only once. These installations can be shared across the team.


-----------------------------------------------------
Requirements, Dependencies, and Supported Compilers
-----------------------------------------------------

List basic requirements, such as language standards (C++, Fortran), CMake
version, etc.

List the compilers we support; i.e., those we regularly build and test with...

List any other dependencies folks need to know about...


.. _tplbuild-label:

----------------------------------------------
Building and Installing Third-party Libraries
----------------------------------------------

We use the `Spack Package Manager <https://github.com/scalability-llnl/spack>`_ 
to manage and build TPL dependencies for the Toolkit. To make the TPL process
easier (you don't really need to learn much about Spack) and automatic, we 
drive it with a python script called ``uberenv.py``, which located in the 
directory 'scripts/uberenv'. Running this script does several things:

  * Clone the Spack repo from GitHub (the most recent version we have tested 
    and we know works!)
  * Perform set up tasks for each version of the TPLs that will be built 
    (e.g., for each compiler that will be used).
  * Invokes Spack to build all TPLs versions an generate a *host-config* file,
    that captures all details of the configuration and build, for each.

The figure illustrates what the script does.

.. figure:: Uberenv.jpg

The uberenv script is run from the top-level Toolkit directory like this::

    $ python ./scripts/uberenv/uberenv.py --prefix {install path} --spec spec  [ --mirror {mirror path} ]

The 'install path' specifies the directory where the TPLs will be installed. 
The 'spec' argument refers to Spack's specification syntax. Typically, a Spack
spec (that's fun to say, no?) indicates a compiler and version for a build.
You can see some examples of this in the python scripts we use to build 
TPLs for the Toolkit development team on LC platforms at LLNL located in
the directory 'scripts/uberenv/llnl_install_scripts'. For more details, please
see the `Spack Spec Documentation <http://spack.readthedocs.io/en/latest/basic_usage.html#specs-dependencies>_`. The 'mirror' argument indicates the location 
of a location where Spack will place the downloaded code for the TPLs. When
building more than one installation of the TPLs, using a mirror will tell 
Spack to only download the distribution for each once and use that for all
installations. To setup a mirror for Spack, run the following before running
the uberenv.py script::

    $ spack mirror create -d {directory} --dependencies uberenv-asctoolkit

Here, 'directory' is the location of the mirror.


.. _toolkitbuild-label:

--------------------------------------
Building and Installing the CS Toolkit
--------------------------------------

We noted in the previous section that......


We use a CMake-based system to configure and build our code, called *BLT*
(see **add link to BLT docs** for more information). 

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


Run CMake directly
^^^^^^^^^^^^^^^^^^^

You can configure the code by running CMake directly and passing it the 
appropriate arguments. For example, to configure, build and install a release 
build with the gcc compiler, you could pass a host-config file CMake::

   $ mkdir build-gnu-release
   $ cd build-gnu-release
   $ cmake -C ./host-configd/surface-chaos_5_x86_64_ib-gcc@4.9.3.cmake \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gnu-release \
     ../src/
   $ make
   $ make install

Alternatively, you could forego the host-config file entirely and pass all the 
argeuments you need to CMakel; for example:: 

   $ mkdir build-gnu-release
   $ cd build-gnu-release
   $ cmake -DCMAKE_C_COMPILER=/usr/apps/gnu/4.9.3/bin/gcc \
     -DCMAKE_CXX_COMPILER=/usr/apps/gnu/4.9.3/bin/g++ \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gnu-release \
     ... \
     ../src/
   $ make
   $ make install

.. note :: The locations of all required third-party libraries must be 
           provided here. These are encoded in our host-config files.

CMake options
^^^^^^^^^^^^^^^

Need to describe CMake options that users would want to provide....Is this 
table correct and up-to-date?

+-----------------------------------+-------------------------------+--------+
|OPTION                             |   Description                 | Default|
+===================================+===============================+========+
|ENABLE_SHARED_LIBS                 |Build shared libraries.        |        |
|                                   |Default is Static libraries    |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_TESTS                       |Builds unit tests              |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_BOOST                       |Enable Boost                   |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_CODECOV                     |Enable code coverage via gcov  |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_CXX11                       |Enables C++11 language support |  ON    | 
+-----------------------------------+-------------------------------+--------+
|ENABLE_FORTRAN                     |Enables Fortran compiler       |  ON    |
|                                   |support                        |        |
+-----------------------------------+-------------------------------+--------+
|ENABLE_ALL_WARNINGS                |Enable extra compiler warnings |        | 
|                                   |in all build targets           |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_WARNINGS_AS_ERRORS          |Compiler warnings treated as   |        |
|                                   |errors.                        | OFF    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_MPI                         |ENABLE MPI                     | OFF    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_OPENMP                      |ENABLE OpenMP                  | OFF    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_BENCHMARKS                  |ENABLE google benchmark        | OFF    |
+-----------------------------------+-------------------------------+--------+


.. _hostconfig-label:

Host-config files
^^^^^^^^^^^^^^^^^^^

We use *host-config* files to track build configurations we support and 
maintain reproducibility. We maintain a collection of such files in the 
'host-configs' directory for platforms and compilers we support. 
When passed to CMake, using the '-C' option, a host-config file initializes 
the CMake cache with the configuration specified in the file. 

.. note :: Need to describe how users would go about generating new
           host config files if they need to...



--------------------------
Make targets
--------------------------

Our system provides a variety of make targets to build individual Toolkit 
components, documentation, run tests, examples, etc. After running CMake 
(using either the python helper script or directly), you can see a listing of
all evailable targets by passing 'help' to make; i.e.,::

   $ make help

The name of each target should be sufficiently descriptive to indicate
what the target does. For example, to run all tests and make sure the
Toolkit components are build properly, execute the following command::

   $ make test

.. note :: Add a table listing and describing the most common make targets
           users may want to use (see table above for format).


.. _appbuild-label:

------------------------------------------
Compiling and Linking with an Application
------------------------------------------

Fill this in...
