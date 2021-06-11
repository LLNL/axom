.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Configuration and Building
==========================

This section provides information about configuring and building
the Axom software after you have cloned the repository.
The main steps for using Axom are:

  #. Configure, build, and install third-party libraries (TPLs) on which Axom depends.
  #. Build and install Axom component libraries that you wish to use.
  #. Build and link your application with the Axom installation.

Depending on how your team uses Axom, some of these steps, such as
installing the Axom TPLs and Axom itself, may need to be done
only once. These installations can be shared across the team.


Requirements, Dependencies, and Supported Compilers
---------------------------------------------------

Basic requirements:
~~~~~~~~~~~~~~~~~~~

  * C++ Compiler with C++11 support
  * CMake
  * Fortran Compiler (optional)

Supported Compilers
~~~~~~~~~~~~~~~~~~~

Axom supports a wide variety of compilers.
Please see the ``<axom_src>/scripts/spack/configs/<platform>/compilers.yaml``
for an up to date list of the currently supported and tested compilers for
each platform.

.. _dependencies-label:

External Dependencies
~~~~~~~~~~~~~~~~~~~~~~

Axom's dependencies come in two flavors:

* Library: contain code that axom must link against
* Tool:  executables that we use as part of our development process, e.g. to generate documentation and format our code.

Unless otherwise marked, the dependencies are optional.

Library Dependencies
""""""""""""""""""""

================== ==================================== ======================
  Library          Dependent Components                 Build system variable
================== ==================================== ======================
  `Conduit`_       Required: Sidre                      CONDUIT_DIR
  `c2c`_           Optional                             C2C_DIR
  `HDF5`_          Optional: Sidre                      HDF5_DIR
  `Lua`_           Optional: Inlet                      LUA_DIR
  `MFEM`_          Optional: Quest                      MFEM_DIR
  `RAJA`_          Optional: Mint, Spin, Quest          RAJA_DIR
  `SCR`_           Optional: Sidre                      SCR_DIR
  `Umpire`_        Optional: Core, Spin, Quest          UMPIRE_DIR
================== ==================================== ======================

.. _Conduit: https://llnl-conduit.readthedocs.io/en/latest
.. _c2c: https://rzlc.llnl.gov/c2c
.. _HDF5: https://www.hdfgroup.org/solutions/hdf5/
.. _Lua: https://www.lua.org/
.. _MFEM: https://mfem.org/
.. _RAJA: https://raja.readthedocs.io/en/main/
.. _SCR: https://computation.llnl.gov/projects/scalable-checkpoint-restart-for-mpi
.. _Umpire: https://umpire.readthedocs.io/en/latest/

Each library dependency has a corresponding build system variable
(with the suffix ``_DIR``) to supply the path to the library's installation directory.
For example, ``hdf5`` has a corresponding variable ``HDF5_DIR``.

.. note::
  Optional `c2c` library is currently only available for configurations on LLNL clusters.


Tool Dependencies
"""""""""""""""""

================== ==================================== ======================
  Tool             Purpose                              Build System Variable
================== ==================================== ======================
  `clangformat`_   Code Style Checks                    CLANGFORMAT_EXECUTABLE
  `CppCheck`_      Static C/C++ code analysis           CPPCHECK_EXECUTABLE
  `Doxygen`_       Source Code Docs                     DOXYGEN_EXECUTABLE
  `Lcov`_          Code Coverage Reports                LCOV_EXECUTABLE
  `Shroud`_        Multi-language binding generation    SHROUD_EXECUTABLE
  `Sphinx`_        User Docs                            SPHINX_EXECUTABLE
================== ==================================== ======================

.. _clangformat: https://releases.llvm.org/10.0.0/tools/clang/docs/ClangFormat.html
.. _CppCheck: http://cppcheck.sourceforge.net/
.. _Doxygen: http://www.doxygen.nl/
.. _Lcov: http://ltp.sourceforge.net/coverage/lcov.php
.. _Shroud: https://shroud.readthedocs.io/en/develop/
.. _Sphinx: http://www.sphinx-doc.org/en/master/

Each tool has a corresponding build system variable (with the suffix ``_EXECUTABLE``)
to supply the tool's executable path. For example, ``sphinx`` has a corresponding build
system variable ``SPHINX_EXECUTABLE``.

.. note::
  To get a full list of all dependencies of Axom's dependencies in an ``uberenv``
  build of our TPLs, please go to the TPL root directory and
  run the following spack command ``./spack/bin/spack spec axom``.


.. _tplbuild-label:


Building and Installing Third-party Libraries
---------------------------------------------

We use the `Spack Package Manager <https://github.com/spack/spack>`_
to manage and build TPL dependencies for Axom. The Spack process works on Linux and macOS
systems. Axom does not currently have a tool to automatically build dependencies for
Windows systems.

To make the TPL process easier (you don't really need to learn much about Spack) and
automatic, we drive it with a python script called ``uberenv.py``, which is located in the
``scripts/uberenv`` directory. Running this script does several things:

  * Clones the Spack repo from GitHub and checks out a specific version
    that we have tested.
  * Configures Spack compiler sets, adds custom package build rules and sets any options
    specific to Axom.
  * Invokes Spack to build a complete set of TPLs for each configuration and generates a
    *host-config* file that captures all details of the configuration and build
    dependencies.

The figure illustrates what the script does.

.. figure:: Uberenv.jpg

The uberenv script is run from Axom's top-level directory like this::

    $ python ./scripts/uberenv/uberenv.py --prefix {install path}  \
                                          --spec spec              \
                                        [ --mirror {mirror path} ]


For more details about ``uberenv.py`` and the options it supports,
see the `uberenv docs <https://uberenv.readthedocs.io/en/latest/>`_

You can also see examples of how Spack spec names are passed to ``uberenv.py``
in the python scripts we use to build TPLs for the Axom development team on
LC platforms at LLNL. These scripts are located in the directory
``scripts/llnl_scripts``.


.. _toolkitbuild-label:


Building and Installing Axom
----------------------------

This section provides essential instructions for building the code.

Axom uses `BLT <https://github.com/LLNL/blt>`_, a CMake-based system, to
configure and build the code. There are two ways to configure Axom:

 * Using a helper script ``config-build.py``
 * Directly invoke CMake from the command line.

Either way, we typically pass in many of the configuration options and variables
using platform-specific *host-config* files.


.. _hostconfig-label:

Host-config files
~~~~~~~~~~~~~~~~~

Host-config files help make Axom's configuration process more automatic and
reproducible. A host-config file captures all build configuration
information used for the build such as compiler version and options,
paths to all TPLs, etc. When passed to CMake, a host-config file initializes
the CMake cache with the configuration specified in the file.

We noted in the previous section that the uberenv script generates a
host-config file for each item in the Spack spec list given to it.
These files are generated by spack in the directory where the
TPLs were installed. The name of each file contains information about the
platform and spec.

For more information, see `BLT's host-config documentation <https://llnl-blt.readthedocs.io/en/develop/tutorial/host_configs.html>`_.


Python helper script
~~~~~~~~~~~~~~~~~~~~

The easiest way to configure the code for compilation is to use the
``config-build.py`` python script located in Axom's base directory;
e.g.,::

   $ ./config-build.py -hc {host-config path}

This script requires that you pass it a *host-config* file. The script runs
CMake and passes it the host-config.
See :ref:`hostconfig-label` for more information.

Running the script, as in the example above, will create two directories to
hold the build and install contents for the platform and compiler specified
in the name of the host-config file.

To build the code and install the header files, libraries, and documentation
in the install directory, go into the build directory and run ``make``; e.g.,::

   $ cd {build directory}
   $ make
   $ make install

.. caution :: When building on LC systems, please don't compile on login nodes.

.. tip :: Most make targets can be run in parallel by supplying the '-j' flag
           along with the number of threads to use.
           E.g. ``$ make -j8`` runs make using 8 threads.

The python helper script accepts other arguments that allow you to specify
explicitly the build and install paths and build type. Following CMake
conventions, we support three build types: 'Release', 'RelWithDebInfo', and
'Debug'. To see the script options, run the script without any arguments;
i.e.,::

   $ ./config-build.py

You can also pass extra CMake configuration variables through the script; e.g.,::

   $ ./config-build.py -hc {host-config file name}          \
                       -DBUILD_SHARED_LIBS=ON               \
                       -DENABLE_FORTRAN=OFF

This will configure cmake to build shared libraries and disable fortran
for the generated configuration.


Run CMake directly
~~~~~~~~~~~~~~~~~~

You can also configure the code by running CMake directly and passing it the
appropriate arguments. For example, to configure, build and install a release
build with the gcc compiler, you could pass a host-config file to CMake::

   $ mkdir build-gcc-release
   $ cd build-gcc-release
   $ cmake -C {host config file for gcc compiler}           \
           -DCMAKE_BUILD_TYPE=Release                       \
           -DCMAKE_INSTALL_PREFIX=../install-gcc-release    \
           ../src/
   $ make
   $ make install

Alternatively, you could forego the host-config file entirely and pass all the
arguments you need, including paths to third-party libraries,
directly to CMake; for example::

   $ mkdir build-gcc-release
   $ cd build-gcc-release
   $ cmake -DCMAKE_C_COMPILER={path to gcc compiler}        \
           -DCMAKE_CXX_COMPILER={path to g++ compiler}      \
           -DCMAKE_BUILD_TYPE=Release                       \
           -DCMAKE_INSTALL_PREFIX=../install-gcc-release    \
           -DCONDUIT_DIR={path/to/conduit/install}          \
           {many other args}                                \
           ../src/
   $ make
   $ make install


CMake Configuration Options
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here are the key build system options in Axom:

+------------------------------+---------+--------------------------------+
| OPTION                       | Default | Description                    |
+==============================+=========+================================+
| AXOM_ENABLE_ALL_COMPONENTS   | ON      | Enable all components          |
|                              |         | by default                     |
+------------------------------+---------+--------------------------------+
| AXOM_ENABLE_<FOO>            | ON      | Enables the axom component     |
|                              |         | named 'foo'                    |
|                              |         |                                |
|                              |         | (e.g. AXOM_ENABLE_SIDRE)       |
|                              |         | for the sidre component        |
+------------------------------+---------+--------------------------------+
| AXOM_ENABLE_DOCS             | ON      | Builds documentation           |
+------------------------------+---------+--------------------------------+
| AXOM_ENABLE_EXAMPLES         | ON      | Builds examples                |
+------------------------------+---------+--------------------------------+
| AXOM_ENABLE_TESTS            | ON      | Builds unit tests              |
+------------------------------+---------+--------------------------------+
| AXOM_ENABLE_TOOLS            | ON      | Builds tools                   |
+------------------------------+---------+--------------------------------+
| BUILD_SHARED_LIBS            | OFF     | Build shared libraries.        |
|                              |         | Default is Static libraries    |
+------------------------------+---------+--------------------------------+
| ENABLE_ALL_WARNINGS          | ON      | Enable extra compiler warnings |
|                              |         | in all build targets           |
+------------------------------+---------+--------------------------------+
| ENABLE_BENCHMARKS            | OFF     | Enable google benchmark        |
+------------------------------+---------+--------------------------------+
| ENABLE_CODECOV               | ON      | Enable code coverage via gcov  |
+------------------------------+---------+--------------------------------+
| ENABLE_FORTRAN               | ON      | Enable Fortran compiler        |
|                              |         | support                        |
+------------------------------+---------+--------------------------------+
| ENABLE_MPI                   | OFF     | Enable MPI                     |
+------------------------------+---------+--------------------------------+
| ENABLE_OPENMP                | OFF     | Enable OpenMP                  |
+------------------------------+---------+--------------------------------+
| ENABLE_WARNINGS_AS_ERRORS    | OFF     | Compiler warnings treated as   |
|                              |         | errors.                        |
+------------------------------+---------+--------------------------------+

If ``AXOM_ENABLE_ALL_COMPONENTS`` is OFF, you must explicitly enable the desired
components (other than 'core', which is always enabled).

See `Axom software documentation <../../../index.html>`_
for a list of Axom's components and their dependencies.

.. note :: To configure the version of the C++ standard, you can supply one of the
           following values for **BLT_CXX_STD**:  'c++11' or 'c++14'.
           Axom requires at least 'c++11', the  default value.

See :ref:`dependencies-label` for configuration variables to specify paths
to Axom's dependencies.


Make targets
------------

Our system provides a variety of make targets to build individual Axom
components, documentation, run tests, examples, etc. After running CMake
(using either the python helper script or directly), you can see a listing of
all available targets by passing 'help' to make; i.e.,::

   $ make help

The name of each target should be sufficiently descriptive to indicate
what the target does. For example, to run all tests and make sure the
Axom components are built properly, execute the following command::

   $ make test


.. _appbuild-label:

Compiling and Linking with an Application
-----------------------------------------

Please see :ref:`using_in_your_project` for examples of how to use Axom in your project.
