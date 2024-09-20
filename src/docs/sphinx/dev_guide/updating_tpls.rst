.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _tpls-label:

*********************
Third-party Libraries
*********************

Axom dependencies are grouped into four categories: Git submodules,
built-in Third-party Libraries (TPLs) in the Axom source tree, system-level
TPLs, and other TPL libraries. The following sections describe how to
install and update these dependencies for Axom.

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    and experimentation?
  * Others?

Determinism
-----------

We strive for as close to deterministic behavior in our builds as possible.
By this, we mean that repeated builds should act the same in the following
regards:

* Set of libraries with their options and versions
* Compilers, compiler flags, and versions
* Installed file and directory structure with permissions


===========================================
Build Scripts and Their Configuration Files
===========================================

There are three levels of build scripts or programs that drive TPL builds.
As you move up the levels, and away from Spack, the scripts require less
configuration and even build multiple sets of TPLs and/or Axom configurations
at a time.

Here is a brief description of what the levels are handling and what important
configuration and input files they use, from lowest level to highest.

Level 1: Spack
--------------

Spack is a multi-platform package manager that builds and installs multiple versions
and configurations of software packages. It has recipes on how to build each package
with variants on each package to customize them to your needs.  For example, Axom
has variants for Fortran and MPI, among others.  These recipes handle how to drive
the individual packages build systems, as well as any packages they depend on.
Spack also handles system level packages, so you can describe where they are on your
system instead of building them from scratch.  You will need to describe which compilers
are available on your system as well.

* Platform specific configuration files live under ``scripts/spack/configs/<platform name>``.
  There is one file (``spack.yaml``) per platform that handles the following:

   * ``compilers``: This section contains the compiler specs that describe the location
     and any other required information about that compiler.  For example, compiler or 
     linker flags.
   * ``packages``: This section describes the system level packages.  For example,
     where they are located and what version they are. This file is very important
     due to its ability to drastically reduce the amount of packages that Spack builds.

* Axom specific Spack package files live under ``scripts/spack/packages``. These override
  the package files that live in Spack's repository here ``var/spack/repos/builtin/packages``.
  We try to minimize these but we have had to alter the existing packages to apply fixes before
  pushing them up to Spack proper or alterations to the recipes that are Axom specific.
  This overriding does not happen at the Spack level, but at the next level, Uberenv.
* `Spack's GitHub repo <https://github.com/spack/spack>`_
* `Spack's documentation <https://spack.readthedocs.io/en/latest/>`_

.. note::
   Spack does not stop at the first error.  It attempts to build as many packages
   as possible.  Due to this, finding the actual error can sometimes be hard but looking
   through the log for a large indented section will help.  The error will
   be in that section and also a message with a path to the full log will be printed
   by Spack afterwards. Searching for ``-build-out.txt`` in your output should
   help.

Level 1: Vcpkg
--------------

Vcpkg is an open-source C++ Library Manager for Windows, Linux, and MacOS by Microsoft.
Axom only uses it for our Windows TPL builds.

* Project specific package files live under ``develop/scripts/vcpkg_ports``.  There are
  two different files for each package:

   * ``portfile.cmake``: This file is the recipe on how to build the package. Vcpkg
     has strict rules about how your project is laid out and you can do the conversion
     in this file as well.
   * ``vcpkg.json``: This is the manifest file that describes information about the
     package.  For example, dependencies, license information, and optional features.

* `Vcpkg's GitHub repo <https://github.com/microsoft/vcpkg>`_
* `Vcpkg's documentation <https://github.com/microsoft/vcpkg#table-of-contents>`_

Level 2: Uberenv
----------------

Uberenv simplifies the use of two level 1 package managers, Spack and Vcpkg.
We rely on Uberenv for two major points: reducing multiple commands into one
and adding as much determinism as possible. The basic workflow in Uberenv is
the following:

#. Setup necessary paths and directories like the base directory where the
   package manager will be installed.
#. Clone the package manager to the specific Git commit.
#. Apply patches to package manager. For example, disabling extra config scopes in Spack.
#. Adds our repositories package repository to Spack, so our packages take precedence.
#. Clean previous temporary information from previous runs that may bleed into this run.
#. Optionally create a package source mirror.
#. Install packages via the selected package manager.

* ``.uberenv_config.json``: This file describes project specific configurations,
  such as, where to download the package manager, what git commit to use, and
  the top level package to install.
* `Uberenv's GitHub repo <https://github.com/LLNL/uberenv>`_
* `Uberenv's documentation <https://uberenv.readthedocs.io/en/latest/>`_

.. note::
   Uberenv's warnings and errors are easy to find by searching the output for ``[ERROR:``
   or ``[Warning:``.  Uberenv will stop at the first error.

Level 3: Build Scripts
----------------------

The file ``axom/scripts/spack/specs.json`` contains a list of all specs
required per platform or machine name. These specs automatically handle
platform differences and contain the full list of compilers and package specs
required to build.

The directory ``axom/scripts/llnl_scripts`` contains three "build" scripts that
are designed to handle building suites of TPLs via Uberenv and Spack.

* ``build_tpls.py``: This script starts by building all TPLs listed in the file
  ``specs.json``. It will generate host-config files and copy them to the base
  of the Axom repository. After building all of the TPLs, it will test Axom
  against those built TPLs as well as test the installed ``using-with-cmake``
  example for correctness. This script stops at the first failed TPL build but
  attempts to build all host-configs against the Axom source with a summary at
  the end of which succeeded or failed.
* ``build_src.py``: This script takes the existing host-configs, or the
  specific one you point at, and builds and tests Axom against them. It also
  tests the ``using-with-cmake`` examples.
* ``build_devtools.py``: This script builds and installs the developer tools
  listed in the ``axom/scripts/spack/packages/axomdevtools/package.py`` Spack
  package. It also uses a different set of Spack configs located in the 
  ``scripts/spack/devtools_config`` directory, so that the regular Spack configs
  can reuse previously built developer tools.

.. note::
   Due to the large amount of information printed to the screen over a full build, the build scripts
   redirect most build step output to log files.  They will not only tell you what command is being run,
   i.e., ``[exe: some/command --with-options]``, but it will tell you the log file being written
   to before it redirects the output from the command, i.e., ``[[log file: /path/to/log``.


=============
Updating TPLs
=============

Git submodules
--------------

Currently, Axom uses three external packages that appear in the repo
as Git submodules. These are the following, including the location of the
package in the Axom source tree:

  * `BLT <https://github.com/LLNL/blt.git>`_, which is the CMake-based build
    system we use. It is located in ``axom/src/cmake/blt``.
  * `Axom Data <https://github.com/LLNL/axom_data.git>`_, which is a collection
    of data files used in testing Axom. It is located in ``axom/data``.
  * `Uberenv <https://github.com/LLNL/uberenv.git>`_, which contains Python
    scripts we use to help automate building third-party dependencies for
    development and deployment. It is located in ``axom/scripts/uberenv``.

There is no software installation process for these dependencies in the 
traditional sense. To update one of these packages in Axom, simply go into
its directory in Axom and check out a new version. If a version is intended
to be changed in the Axom repo, make the version change on a branch and 
submit a GitHub pull request as you would do for other software changes.
More info on :ref:`building-axom-label`.

Built-in TPLs
-------------

Axom uses several lightweight, header-only libraries internally, which are
exposed for downstream customers to use if they wish.

  * `CLI11 <https://github.com/CLIUtils/CLI11>`_ is a command line parser
    for C++11 and beyond that provides a rich feature set with a simple and
    intuitive interface.
  * `fmt <https://github.com/fmtlib/fmt>`_ is an open-source formatting
    library providing a fast and safe alternative to C stdio and C++ iostreams.
  * `sol <https://github.com/ThePhD/sol2>`_ is a C++ library binding to Lua.
  * `Sparsehash <https://github.com/sparsehash/sparsehash>`_ contains several
    hash-map implementations.

.. note:: Axom patches all built-in TPLs to be under the ``axom`` namespace.
   This is to prevent symbol collisions with other projects, either our
   dependencies or downstream customers who wish their own versions.  For
   example, ``fmt::format("foo")`` is ``axom::fmt::format("foo")``.

They can be found in the directory: ``axom/src/thirdparty/axom``. The basic 
instructions on how to update a built-in TPL are as follows:

#. Download the new release and override the source that is already there.
   This can often involve removing files no-longer needed but most of the
   current ones are a single header file.

#. Review and apply the existing patch files. More than likely, you will not
   be able to directly apply the patch but it will give you the general idea
   on what needs to be applied. For example, the namespace update mentioned above.

#. Ensure that the build and tests still pass. More info on :ref:`testing-label`.

#. Follow the normal pull request work flow. More info on :ref:`pullrequest-label`.

.. _local-tpls-label:

Local Third-party Library Installation
--------------------------------------

It is often useful to have a different set of TPLs during the development process.
For example, you may want to try out a new library or version of an existing library.

From the top-level Axom directory, run the following script to build all TPLs
for all existing compiler specs on the platform you are currently on::

$ ./scripts/llnl_scripts/build_tpls.py -d local/install/path

where ``local/install/path`` is a directory location where you want the 
libraries to be installed.

It will output whether the TPL install succeeded and, 
subsequently, whether an Axom build against the TPL install succeeded.

Running the script produces new host-config files (i.e., CMake cache files) 
that you can use to build and test Axom with the installation, if issues
arise. The generated host-config files will be located in the top-level Axom
directory of your local clone of the repo. If any changes to Axom code are 
needed to work with the TPL update(s), make the changes and test them.

.. note:: You can build a subset of TPLs for a platform, by using
          the ``uberenv.py`` script in the top-level Axom directory.
          For example:: 

            python3 ./scripts/uberenv/uberenv.py --prefix /my/tpl/path --spec clang@10.0.0~cpp14+devtools+mfem+c2c

          will build the TPLs for the clang 10.0.0 compiler, install them
          to the ``/my/tpl/path`` directory, and generate a host-config file
          that you can use to build Axom and its tests. Please see the
          ``scripts/spack/specs.json`` file for a current list of tested specs. 


Shared Third-party Library Installation Steps
---------------------------------------------

The following instructions describe how to install local copies of Axom
TPLs on Livermore Computing (LC) platforms and recreate our Docker containers
with a new set of TPLs. Typically, this process is followed when you want to 
update one or more TPLs which Axom depends on. After they are built and
the required changes are merged into develop, they will be available for
other Axom developers to use during development, in Axom GitLab CI testing, etc.

#. **Working on a local branch.** 
   Make a local clone of the Axom repo and create a branch to work on.

#. **Changing versions of system packages or other TPLs.**
   To change a version of a system package, which applies to an LC platforms 
   or a Docker container image we use for CI testing on GitHub, go into
   the directory ``axom/scripts/spack/configs``. There you will find a 
   sub-directory for each supported LC system type. Each sub-directory
   has a ``spack.yaml`` file which contains an entry for each system level
   package we rely on. Find the entry for the library you wish to update and 
   change the version number. Do this for each system you want to test/change,
   including configurations in the ``docker`` subdirectory.

   .. note:: Inside of the ``spack.yaml`` for each system package directory,
             there is a ``compilers`` section containing compiler and 
             version information for compilers we use for development and 
             testing. If you wish to test and build with a new compiler or 
             version on a system, modify the appropriate ``spack.yaml`` 
             file.

   To change a version of a non-system TPL, go into the 
   ``axom/scripts/spack/configs`` directory. There you will find a 
   sub-directory for each system we test on which contains a Spack
   package file ``package.py``. TPL versions are pinned in those package files.
   Modify the contents of the Spack package file ``package.py`` in each
   package sub-directory as needed to change TPL version numbers.

   .. note:: Before continuing, you should test that the installation works
             on all LC systems with the steps in :ref:`local-tpls-label`.


#. **Install TPLs on all required LC machines.**
   This step needs to be run on each of the machines named in Axom's standard host-configs.
   When you are confident that everything is correct, become the service user
   ``atk`` via the following command::

   $ xsu atk

   .. note:: This command requires special access permissions. If you need them, contact the Axom team.

   Run the corresponding command for the system you are on::

     # blueos
     $ lalloc 1 -W 120 scripts/llnl_scripts/build_tpls.py
     
     # toss_4
     $ srun -N1 --interactive -t 120 scripts/llnl_scripts/build_tpls.py

   This script will build all third-party libraries for all compilers specs
   for the machine you are on. These will be installed into the shared LC directory
   ``/usr/workspace/axom/libs/<SYS_TYPE>/<time date>/<compiler>``
   used by Axom developers. When completed, they will produce new host-config
   files for each configuration. These host-configs will be at the base of the repository
   and named in the following pattern: ``<machine name>-<SYS_TYPE>-<compiler spec>.cmake``
   Give these files to your regular user account
   and log back in to that account. Copy these new host-config files to the
   ``host-configs`` subdirectory and commit them to your branch. Make sure all
   file changes from all previous steps are also committed and pushed upstream.

   .. note:: If this step fails, delete the time date stamped directory that was created.
             If you forget to do this, it will eventually be deleted by hand in bulk when
             they are past a certain age and no longer needed.

#. **Build new Docker images.**
   We utilize Docker images that have pre-built TPLs in our GitHub CI checks.
   To build these, go to our
   `GitHub Actions <https://github.com/LLNL/axom/actions/workflows/docker_build_tpls.yml>`_
   page. Click on "Actions" and then on "Docker TPL build" in the "Workflows" menu.
   Find the "Run Workflow" drop-down menu, select your branch, and click on the "Run workflow"
   button. This will launch the build of the docker images.

   When the docker image build completes, click on your build and find the
   "Artifacts" listed at the bottom of the page. These contain host-configs
   for building Axom on the docker images. Download them and copy them to
   Axom's ``host-configs/docker`` subdirectory. Rename them to match the corresponding
   host-config.

#. **Update Azure Pipelines to the new Docker images.**
   To complete the setup of the new docker images, the ``Compiler_ImageName``
   entries in ``azure-pipelines.yaml`` at the top-level directory must be updated
   with the timestamped names of the new images. The new names can be found in
   the log files from the successful GitHub action. On the left of the page for
   the successful action is a "Jobs" menu. Click on each job and then find
   the "Get dockerhub repo name" section of the log. The second line of the
   section there should be an entry of the form ``axom/tpls:clang-10_12-18-20_00h-10m``.
   Copy the name beginning with ``axom/tpls`` to the appropriate locations
   in ``azure-pipelines.yaml``. Repeat this with the names from each compiler
   job used in the GitHub action. 
   Axom's docker images are hosted on our `DockerHub <https://hub.docker.com/r/axom/tpls/tags>`_ page.

#. Make sure all changes in your branch are committed and pushed, and create
   a pull request for a merge to develop. If everything went well, all checks
   on your GitHub PR should pass.
 
