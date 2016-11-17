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

Basic requirements:

  * C++ 98/11
  * CMake 3.1.2
  * Fortran (Optional)


Compilers we support:

  * Clang - 3.5.0
  * GCC - 4.7.1, 4.9.3
  * IBM BGQOS - 12.1.012a
  * Intel - 15.0.187, 16.0.109
  

Package Dependencies:

.. danger::
  Instead of this exhaustive list (which will go stale fast), the more important things to list are the 
  the TPLs we directly use and the  **ZZZ_DIR**, **ZZZ_EXECUTABLE**, etc  
  cmake args we use to pull those TPLS in our build system. For example, folks really don't need to know that Bison is used, we have no logic for it in our build system. It actually it may only be needed due to some spack dependency chain that could change over time. If people do want to know this info, we could point them 
  a few spack commands that could show whatever the current deps are (I don't know what those commands are off hand)


  * Bison 3.0.4
  * Boost 1.58.0
  * bzip2 1.0.6
  * Conduit 2016-05-18
  * Doxygen 1.8.11
  * Flex 2.6.0
  * HDF5 1.8.16
  * lcov-1.11
  * libsigsegv 2.10
  * lua 5.1.5
  * m4 1.4.17
  * ncurses 6.0
  * openssl 1.0.2h
  * py-alabaster 0.7.7
  * py-babel 1.3
  * py-breathe 4.0.0
  * py-cogapp 2.4
  * py-docutils 0.12
  * py-jinja2 2.8
  * py-markupsafe 0.23
  * py-parsley 1.2
  * py-pygments 2.1
  * py-pyyaml 3.11
  * py-setuptools 18.1
  * py-six 1.9.0
  * py-snowballstemmer 1.2.1
  * py-sphinx-rtd-theme 0.1.9
  * py-sphinx 1.3.6
  * py-tz 2015.7
  * python 2.7.11
  * readline 6.3
  * sparsehash-headers 2.0.2
  * sqlite 3.8.5
  * uberenv-asctoolkit 0.1
  * uncrustify 0.61
  * zlib 1.2.8


.. danger::
  Here is one possible format, we can also add urls and version info



================== ==================================== 
  Library            Dependent Components
================== ==================================== 
  HDF5               Sidre, Spio
  Conduit            Sidre, Spio
  Sparse Hash        Sidre (optional)
  Boost              Slam, Quest
================== ==================================== 

================== ==================================== 
  Tool               Purpose
================== ==================================== 
  Shroud            Multi-language binding generation 
  Doxygen            Source Code Docs
  Sphinx             User Docs
  Breathe            Doxygen integration for Sphinx
  Uncrustify         Code Style Checks 
  Lcov               Code Coverage Reports
================== ==================================== 



.. _tplbuild-label:

----------------------------------------------
Building and Installing Third-party Libraries
----------------------------------------------

We use the `Spack Package Manager <https://github.com/scalability-llnl/spack>`_ 
to manage and build TPL dependencies for the Toolkit. To make the TPL process
easier (you don't really need to learn much about Spack) and automatic, we 
drive it with a python script called ``uberenv.py``, which located in the 
directory ``scripts/uberenv``. Running this script does several things:

  * Clones the Spack repo from GitHub and checks out a specific version that we have tested.
  * Configures Spack compiler sets, adds custom package build rules and set any options specific to the Toolkit. 
  * Invokes Spack to build a complete set of TPLs for a configuration an generates a *host-config* file,
    that captures all details of the configuration and the built dependencies.

The figure illustrates what the script does.

.. figure:: Uberenv.jpg

The uberenv script is run from the top-level Toolkit directory like this::

    $ python ./scripts/uberenv/uberenv.py --prefix {install path} --spec spec  [ --mirror {mirror path} ]

Here is a break down of the options that control how ``uberenv.py`` builds dependencies:

 ================== ==================================== ======================================
  Option             Description                          Default
 ================== ==================================== ======================================
  --prefix           Destination directory                ``uberenv_libs``
  --spec             Spack spec                           linux: **%gcc**
                                                          osx: **%clang**
  --compilers-yaml   Spack compilers settings file        ``scripts/uberenv/compilers.yaml``
  --mirror           Spack source mirror location         **Not used**
  --create-mirror    Establish a new Spack source mirror  ``False``
                     with out building TPLs
 ================== ==================================== ======================================

Default invocation on Linux:

.. code:: bash

    python scripts/uberenv/uberenv.py --prefix uberenv_libs \
                                      --spec %gcc \
                                      --compilers-yaml scripts/uberenv/compilers.yaml

Default invocation on OSX:

.. code:: bash

    python scripts/uberenv/uberenv.py --prefix uberenv_libs \
                                      --spec %clang \
                                      --compilers-yaml scripts/uberenv/compilers.yaml
 

The 'install path' specifies the directory where the TPLs will be installed. 
The 'spec' argument refers to Spack's specification syntax. Typically, a Spack
spec (that's fun to say, no?) indicates the specific version of a specific compiler to use for the build.
We manage the set of compilers Spack supports in the ``scripts/uberenv/compilers.yaml`` file. 

You can edit ``scripts/uberenv/compilers.yaml`` or use the **--compilers-yaml** option to select another file to set the  compiler settings used by Spack. See the `Spack Compiler Configuration <http://spack.readthedocs.io/en/latest/getting_started.html#manual-compiler-configuration>`_
documentation for details.

For OSX, the defaults in ``compilers.yaml`` are X-Code's clang and gfortran from https://gcc.gnu.org/wiki/GFortranBinaries#MacOS. 

.. note::
    uberenv.py forces Spack to ignore ``~/.spack/compilers.yaml`` to avoid conflicts
    and surprises from a user's specific Spack settings on HPC platforms.


You can also see examples of how Spack spec names are passed to ``uberenv.py`` in the python scripts we use to build 
TPLs for the Toolkit development team on LLNL's LC platforms. These scripts are located in
the directory ``scripts/uberenv/llnl_install_scripts``. 

The 'mirror' argument provides a location for Spack to store the downloaded source code for TPL dependencies. When
building more than one installation of the TPLs, using a mirror will allow Spack to skip downloads for source code that was already already obtained during a prior build. 

When the 'create-mirror' argument is used, ``uberenv.py`` establishes a Spack mirror and downloads the source for all TPL dependencies into this mirror. It does not build any TPLs. This option is used to obtain a copy of source code for all necessary TPLs so it can be transfered to another system for builds.

.. _toolkitbuild-label:

--------------------------------------
Building and Installing the CS Toolkit
--------------------------------------

We use a CMake-based system, *BLT*, to configure and build the Toolkit
(see **add link to BLT docs** for more information). This section 
provides essential instructions for building the code.


.. _hostconfig-label:

Host-config files
^^^^^^^^^^^^^^^^^^^

We use host-config files to make building the Toolkit more automatic and
easily reproducible. A host-config file captures all build configuration 
information used for the build such as compiler version and options, 
paths to all TPLs, etc. When passed to CMake, a host-config file initializes
the CMake cache with the configuration specified in the file.

We noted in the previous section that the uberenv script generates a
'host-config' file for each item in the Spack spec list given to it.
These files are located in the directory ``spack/bin/spack`` where the
TPLs were installed. The name of each file contains information about the
platform and spec.


Python helper script
^^^^^^^^^^^^^^^^^^^^^

The easiest way to configure the code for compilation is to use the 
``config-build.py`` python script in the 'scripts' directory; 
e.g.,::

   $ ./scripts/config-build.py -hc {host-config file name}

This script requires that you pass it a *host-config* file. The script runs 
CMake and passes it the host-config. See :ref:`hostconfig-label` 
for more information.

Running the script, as in the example above, will create two directories to 
hold the build and install contents for the platform and compiler specified 
in the name of the host-config file. 

To build the code and intall the header files, libraries, and documentation 
in the install directory, go into the build directory and run ``make``; e.g.,::

   $ cd {build directory}
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

You can also configure the code by running CMake directly and passing it the 
appropriate arguments. For example, to configure, build and install a release 
build with the gcc compiler, you could pass a host-config file to CMake::

   $ mkdir build-gcc-release
   $ cd build-gcc-release
   $ cmake -C {host config file for gcc compiler} \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gcc-release \
     ../src/
   $ make
   $ make install

Alternatively, you could forego the host-config file entirely and pass all the 
argeuments you need to CMake; for example:: 

   $ mkdir build-gcc-release
   $ cd build-gcc-release
   $ cmake -DCMAKE_C_COMPILER={path to gcc compiler} \
     -DCMAKE_CXX_COMPILER={path to g++ compiler} \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_INSTALL_PREFIX=../install-gcc-release \
     {many other args} \
     ../src/
   $ make
   $ make install

.. note :: The locations of all required third-party libraries must be 
           provided here. These are encoded in our host-config files.

CMake options
^^^^^^^^^^^^^^^

.. note :: Summarize (in table) CMake options that users may want to provide
           Check what's there now for correctness.

+-----------------------------------+-------------------------------+--------+
|OPTION                             |   Description                 | Default|
+===================================+===============================+========+
|ENABLE_ALL_COMPONENTS              |Enables all components         |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_ALL_WARNINGS                |Enable extra compiler warnings |  ON    | 
|                                   |in all build targets           |        |
+-----------------------------------+-------------------------------+--------+
|ENABLE_BENCHMARKS                  |Enable google benchmark        |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_BOOST                       |Enable Boost                   |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_CFORTRAN_API                |Enable C to Fortran interface  |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_CODECOV                     |Enable code coverage via gcov  |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_CXX11                       |Enable C++11 language support  |  ON    | 
+-----------------------------------+-------------------------------+--------+
|ENABLE_FORTRAN                     |Enable Fortran compiler        |  ON    |
|                                   |support                        |        |
+-----------------------------------+-------------------------------+--------+
|ENABLE_MPI                         |Enable MPI                     |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_OPENMP                      |Enable OpenMP                  |  OFF   |
+-----------------------------------+-------------------------------+--------+
|ENABLE_SHARED_LIBS                 |Build shared libraries.        |  OFF   |
|                                   |Default is Static libraries    |        |
+-----------------------------------+-------------------------------+--------+
|ENABLE_TESTS                       |Builds unit tests              |  ON    |
+-----------------------------------+-------------------------------+--------+
|ENABLE_WARNINGS_AS_ERRORS          |Compiler warnings treated as   |  OFF   |
|                                   |errors.                        |        |
+-----------------------------------+-------------------------------+--------+

.. danger:: 
    We are going to pull out boost, but if we kept it in - we should simply support
    via BOOST_DIR, not ENABLE_BOOST + BOOST_ROOT
    


CMake Options used to include Third-party Libraries:

+-----------------------------------+---------------------------------------------------+
|OPTION                             |   Description                                     |
+===================================+===================================================+
|HDF5_DIR                           | Path to HDF5 install                              |
+-----------------------------------+---------------------------------------------------+
|CONDUIT_DIR                        | Path to Conduit install                           |
+-----------------------------------+---------------------------------------------------+
|PYTHON_EXECUTABLE                  | Path to Python executable                         |
+-----------------------------------+---------------------------------------------------+
|SPARSEHASH_DIR                     | Path to Sparsehash install                        |
+-----------------------------------+---------------------------------------------------+

CMake Options used to enable Software Development Tools (should these go in BLT docs and link here?):

+-----------------------------------+---------------------------------------------------+
|OPTION                             |   Description                                     |
+===================================+===================================================+
|SPHINX_EXECUTABLE                  | Path to sphinx-build executable (support via BLT) |
+-----------------------------------+---------------------------------------------------+
|DOXYGEN_EXECUTABLE                 | Path to doxygen executable (support via BLT)      |
+-----------------------------------+---------------------------------------------------+
|UNCRUSTIFY_EXECUTABLE              | Path to uncrustify executable (support via BLT)   |
+-----------------------------------+---------------------------------------------------+


.. danger::
    TODO: LCOV_PATH, GENHTML_PATH, GCOV_PATH  -- these aren't named consistently (_EXECUTABLE suffix?) 



--------------------------
Make targets
--------------------------

Our system provides a variety of make targets to build individual Toolkit 
components, documentation, run tests, examples, etc. After running CMake 
(using either the python helper script or directly), you can see a listing of
all available targets by passing 'help' to make; i.e.,::

   $ make help

The name of each target should be sufficiently descriptive to indicate
what the target does. For example, to run all tests and make sure the
Toolkit components are built properly, execute the following command::

   $ make test

.. note :: Add a table listing and describing the most common make targets
           users may want to use (see table above for format).


.. _appbuild-label:

------------------------------------------
Compiling and Linking with an Application
------------------------------------------

Fill this in...
