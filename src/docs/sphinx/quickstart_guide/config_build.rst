.. ##
.. ## Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

==========================
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


---------------------------------------------------
Requirements, Dependencies, and Supported Compilers
---------------------------------------------------

.. danger:: This needs to be updated and organized better...

Basic requirements:

  * C++ 98/11
  * CMake 3.1.2
  * Fortran (optional)

Compilers we support (listed with minimum supported version):

  * Clang 3.9.1
  * GCC 4.9.3
  * IBM XL 13
  * Intel 15
  * Microsoft Visual Studio 2015
  * Microsoft Visual Studio 2015 with the Intel toolchain

Please see the ``<axom>/scripts/uberenv/compilers.yaml`` for an up to date
list of the supported compiler version numbers on each platform. 

External Dependencies:

Axom has the following external dependencies. 
Unless otherwise marked, the dependencies are optional.
  
  * conduit 0.3.1 [Required for Sidre]
  * doxygen 1.8.11
  * hdf5 1.8.16
  * lcov 1.11
  * mfem 3.3.2
  * python 2.7.15
  * py-shroud 0.10.0
  * py-sphinx 1.3.6
  * scr 1.2.1
  * uncrustify 0.61

Each dependency in the above list has a corresponding variable that can be 
supplied to the build system. These variables either list a path to the
installation directory  (in which case the variable has the suffix ``_DIR``)
or has the path to an executable (with the ``_EXECUTABLE`` suffix).
For example, ``hdf5`` has a corresponding variable ``HDF5_DIR``
and ``python`` has a corresponding build system variable ``PYTHON_EXECUTABLE``.

.. note::
  To get a full list of all dependencies of Axom's dependencies in an ``uberenv``
  build of our TPLs, please go to the TPL root directory and 
  run the following spack command ``./spack/bin/spack find``.

.. danger::
  Here is one possible format, we can also add urls and version info


================== ====================================
  Library            Dependent Components
================== ====================================
  HDF5               Sidre
  Conduit            Sidre
  SCR                Sidre (optional)
================== ====================================

================== ====================================
  Tool               Purpose
================== ====================================
  Shroud             Multi-language binding generation
  Doxygen            Source Code Docs
  Sphinx             User Docs
  Breathe            Doxygen integration for Sphinx
  Uncrustify         Code Style Checks
  Lcov               Code Coverage Reports
================== ====================================



.. _tplbuild-label:

---------------------------------------------
Building and Installing Third-party Libraries
---------------------------------------------

We use the `Spack Package Manager <https://github.com/scalability-llnl/spack>`_
to manage and build TPL dependencies for Axom. To make the TPL process
easier (you don't really need to learn much about Spack) and automatic, we
drive it with a python script called ``uberenv.py``, which is located in the
``scripts/uberenv`` directory. Running this script does several things:

  * Clones the Spack repo from GitHub and checks out a specific version that we have tested.
  * Configures Spack compiler sets, adds custom package build rules and sets any options specific to Axom.
  * Invokes Spack to build a complete set of TPLs for each configuration and generates a *host-config* file that captures all details of the configuration and build dependencies.

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
spec ("Spack spec" that's fun to say, huh?) indicates a specific version of
a particular compiler to use for the build. We manage the set of compilers
we support in the ``scripts/uberenv/compilers.yaml`` file.

You can edit ``scripts/uberenv/compilers.yaml`` or use the **--compilers-yaml**
option to select another file to set the  compilers and setting you want.
See the `Spack Compiler Configuration <http://spack.readthedocs.io/en/latest/getting_started.html#manual-compiler-configuration>`_ documentation for details.

For OSX, the defaults in ``compilers.yaml`` are X-Code's clang and gfortran
from `X-code for OSX <https://gcc.gnu.org/wiki/GFortranBinaries#MacOS>`_.

.. note::
    uberenv.py forces Spack to ignore ``~/.spack/compilers.yaml`` to avoid
    conflicts and surprises from a user's specific Spack settings.


You can also see examples of how Spack spec names are passed to ``uberenv.py``
in the python scripts we use to build TPLs for the Axom development team on
LC platforms at LLNL. These scripts are located in the directory
``scripts/uberenv/llnl_install_scripts``.

The 'mirror' argument provides a location for Spack to store the downloaded
source code for TPL dependencies. When building more than one installation
of the TPLs, using a mirror will allow Spack to skip downloads for source
code that was already obtained during a prior build.

When the 'create-mirror' argument is used, ``uberenv.py`` establishes a Spack
mirror and downloads the source for all TPL dependencies into this mirror.
It does not build any TPLs. This option is used to obtain a copy of source
code for all necessary TPLs so it can be transferred to another system for
building there.


.. _toolkitbuild-label:

----------------------------
Building and Installing Axom
----------------------------

We use a CMake-based system, called `BLT <https://github.com/LLNL/blt>`_, to
configure and build Axom. This section provides essential instructions for
building the code.

.. note:: Add instructions for "developer" builds vs. "user" builds.


.. _hostconfig-label:

Host-config files
^^^^^^^^^^^^^^^^^

We use host-config files to make building Axom more automatic and
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
^^^^^^^^^^^^^^^^^^^^

The easiest way to configure the code for compilation is to use the
``config-build.py`` python script in the base directory;
e.g.,::

   $ ./config-build.py -hc {host-config file name}

This script requires that you pass it a *host-config* file. The script runs
CMake and passes it the host-config. See :ref:`hostconfig-label`
for more information.

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

   $ ./config-build.py -hc {host-config file name} \
                       -DAXOM_ENABLE_PYTHON=ON -DENABLE_FORTRAN=OFF

This will enable python and disable fortran for the generated configuration.


Run CMake directly
^^^^^^^^^^^^^^^^^^

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
arguments you need directly to CMake; for example::

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
^^^^^^^^^^^^^

.. note :: Summarize (in table) CMake options that users may want to provide
           Check what's there now for correctness. Move options for developers
           into separate table here (for convenience) or to Dev Guide?

+------------------------------+--------------------------------+---------+
| OPTION                       | Description                    | Default |
+==============================+================================+=========+
| AXOM_ENABLE_ALL_COMPONENTS   | Enable all components          | ON      |
|                              | by default                     |         |
+------------------------------+--------------------------------+---------+
| ENABLE_ALL_WARNINGS          | Enable extra compiler warnings | ON      |
|                              | in all build targets           |         |
+------------------------------+--------------------------------+---------+
| ENABLE_BENCHMARKS            | Enable google benchmark        | OFF     |
+------------------------------+--------------------------------+---------+
| ENABLE_CODECOV               | Enable code coverage via gcov  | ON      |
+------------------------------+--------------------------------+---------+
| ENABLE_FORTRAN               | Enable Fortran compiler        | ON      |
|                              | support                        |         |
+------------------------------+--------------------------------+---------+
| ENABLE_MPI                   | Enable MPI                     | OFF     |
+------------------------------+--------------------------------+---------+
| ENABLE_OPENMP                | Enable OpenMP                  | OFF     |
+------------------------------+--------------------------------+---------+
| ENABLE_SHARED_LIBS           | Build shared libraries.        | OFF     |
|                              | Default is Static libraries    |         |
+------------------------------+--------------------------------+---------+
| AXOM_ENABLE_TESTS            | Builds unit tests              | ON      |
+------------------------------+--------------------------------+---------+
| AXOM_ENABLE_DOCS             | Builds documentation           | ON      |
+------------------------------+--------------------------------+---------+
| AXOM_ENABLE_EXAMPLES         | Builds examples                | ON      |
+------------------------------+--------------------------------+---------+
| ENABLE_WARNINGS_AS_ERRORS    | Compiler warnings treated as   | OFF     |
|                              | errors.                        |         |
+------------------------------+--------------------------------+---------+

If 'AXOM_ENABLE_ALL_COMPONENTS' is OFF, you must explicitly enable the desired
components (other than 'common', which is always enabled).

.. note :: To configure the version of the C++ standard, you can supply one of the
           following values for **BLT_CXX_STD**:  'c++98', 'c++11' or 'c++14'.
           The default is 'c++11'.


CMake Options used to include Third-party Libraries:

+-------------------+-------------------------------+
| OPTION            | Description                   |
+===================+===============================+
| HDF5_DIR          | Path to HDF5 install          |
+-------------------+-------------------------------+
| CONDUIT_DIR       | Path to Conduit install       |
+-------------------+-------------------------------+
| MFEM_DIR          | Path to MFEM install          |
+-------------------+-------------------------------+
| PYTHON_EXECUTABLE | Path to Python executable     |
+-------------------+-------------------------------+


CMake Options used to enable Software Development Tools (should these go in BLT docs and link here?):

+-----------------------+---------------------------------------------------+
| OPTION                | Description                                       |
+=======================+===================================================+
| SPHINX_EXECUTABLE     | Path to sphinx-build executable (support via BLT) |
+-----------------------+---------------------------------------------------+
| DOXYGEN_EXECUTABLE    | Path to doxygen executable (support via BLT)      |
+-----------------------+---------------------------------------------------+
| UNCRUSTIFY_EXECUTABLE | Path to uncrustify executable (support via BLT)   |
+-----------------------+---------------------------------------------------+


.. danger::
    TODO: LCOV_PATH, GENHTML_PATH, GCOV_PATH  -- aren't named consistently (_EXECUTABLE suffix?)



------------
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

.. note :: Add a table listing and describing the most common make targets
           users may want to use (see table above for format).


.. _appbuild-label:

-----------------------------------------
Compiling and Linking with an Application
-----------------------------------------

Incorporating Axom as a Git-Submodule to a CMake-Based Application
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you are working on a project based on CMake_
you may want to incorporate Axom as Git submodule as follows:

1. Add Axom as a git submodule to your project, for example: ::

   $ git submodule add ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git <path/to/axom>

.. note::
      If you are not using BLT_ in your project, you'll have to issue the
      following: ::

         git submodule update --init --recursive

      This will put BLT_ in `axom/src/cmake/blt`.

2. Add the following line in the associated "CMakeLists.txt" for your project: ::

      add_subdirectory( axom )

.. _CMake: https://cmake.org
.. _BLT: https://github.com/LLNL/blt
