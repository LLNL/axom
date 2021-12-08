.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _tpls-label:

*********************
Third-party Libraries
*********************

Axom dependencies are grouped into four categories: Git submodules,
built-in TPLs in the Axom source tree, system-level TPLs, and other 
TPL libraries. The following sections describe how to install and update 
these dependencies for Axom.

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    and experimentation?
  * Others?


=============
Updating TPLs
=============

Git submodules
--------------

Currently, Axom uses three external packages that appear in the repo
as Git submodules. These are the following, including the location of the
package in the Axom source tree:

  * `BLT <https://github.com/LLNL/blt.git>`_, which is the CMake-based build
    system we use. Location: ``axom/src/cmake/blt``.
  * `Axom Data <https://github.com/LLNL/axom_data.git>`_, which is where we
    maintain data files used in testing Axom. Location: ``axom/data``.
  * `Uberenv <https://github.com/LLNL/uberenv.git>`_, which contains Python
    scripts we use to help automate building third-party dependencies for
    development and deployment. Location: ``axom/scripts/uberenv``.

There is no software installation process for these dependencies in the 
traditional sense. To update one of these packages in Axom, simply go into
its directory in Axom and check out a new version. If a version is intended
to be changed in the Axom repo, make the version change on a branch and 
submit a GitHub pull request as you would do for other software changes.
More info on :ref:`toolkitbuild-label`.

Built-in TPLs
-------------

Axom several lightweight header-only libraries that we use internally and
expose for downstream customers to use if they wish.

  * `CLI11 <https://github.com/CLIUtils/CLI11>`_, is a command line parser
    for C++11 and beyond that provides a rich feature set with a simple and
    intuitive interface.
  * `fmt <https://github.com/fmtlib/fmt>`_, is an open-source formatting
    library providing a fast and safe alternative to C stdio and C++ iostreams.
  * `sol <https://github.com/ThePhD/sol2>`_, is a C++ library binding to Lua.
  * `Sparsehash <https://github.com/sparsehash/sparsehash>`_, contains several
    hash-map implementations.

.. note:: Axom patches all built-in TPLs to be under the ``axom`` namespace.
   This is to prevent symbol collisions with other projects, either our
   dependencies or downstream customers who wish their own versions.  For
   example, ``fmt::format("foo")`` is ``axom::fmt::format("foo")``.

They can be found in the directory: ``src/thirdparty``. The basic 
instructions on how to update a built-in TPL are as follows:

#. Download the new release and override the source that is already there.
   This can often involve removing files no-longer needed but most of the
   current ones are a single header file.

#. Review and apply the existing patch files. More than likely, you will not
   be able to directly apply the patch but it will give you the general idea
   on what needs to be applied.  For example, the namespace update mentioned above.

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

.. note:: You can build a subset of TPLs for a platform, by passing a Spack
          spec arguments to the `build_tpls.py` script. For example,

          ``--spec clang@10.0.0~cpp14+devtools+mfem+c2c``

          will build the TPLs for the clang 10.0.0 compiler. Please see the
          ``scripts/spack/specs.json`` file for a list of currently tested specs. 


Shared Third-party Library Installation Steps
---------------------------------------------

The following instructions describe how to install local copies of Axom
TPLs on Livermore Computing (LC) platforms and recreate our Docker containers
with a new set of TPLs. Typically, this process is followed when you want to 
update one or more TPLs which Axom depends on. After they are built and
the required changes are merged into develop, they will be available for
other Axom developers to use during development, in Axom Gitlab CI testing, etc.

#. **Working on a local branch.** 
   Make a local clone of the Axom repo and create a branch to work on.

#. **Changing versions of system packages or other TPLs.**
   To change a version of a system package, which applies to an LC platforms 
   or a Docker container image we use for CI testing on GitHub, go into
   the directory ``axom/scripts/spack/configs``. There you will find a 
   sub-directory for each supported LC system type. Each sub-directory
   has a ``packages.yaml`` file which contains an entry for each system level
   package we rely on. Find the entry for the library you wish to update and 
   change the version number. Do this for each system you want to test/change,
   including configurations in the ``docker`` subdirectory.

   .. note:: Alongside each ``packages.yaml`` in each system package directory,
             there is a ``compilers.yaml`` file containing compiler and 
             version information for compilers we use for development and 
             testing. If you wish to test and build with a new compiler or 
             version on a system, modify the appropriate ``compilers.yaml`` 
             file.

   To change a version of a non-system TPL, go into the 
   ``axom/scripts/spack/packages`` directory. There you will find a 
   sub-directory for each TPL Axom uses. Modify the contents of the Spack
   package file ``package.py`` in each package sub-directory as needed. 

   .. note:: Before continuing, you should test that the installation works
             on all LC systems with the steps in :ref:`local-tpls-label`.


#. **Install TPLs on all required LC machines.**
   This step needs to be run on each of the machines named in Axom's standard host-configs.
   When you are confident that everything is correct, become the service user
   ``atk`` via the following command::

   $ xsu atk

   .. note:: This command requires a certain level of permission.

   Run the corresponding command for the system you are on::

     # blueos
     $ lalloc 1 -W 120 scripts/llnl/build_tpl.py
     
     # toss3
     $ srun -N1 --interactive -t 120 scripts/llnl/build_tpl.py

   This script will build all third-party libraries for all compilers specs
   for the machine you are on. These will be installed into shared directories
   used by Axom developers. When completed, they will produce new host-config
   files for each configuration. These host-configs will be at the base of the repository
   and named in the following pattern: ``<machine name>-<SYS_TYPE>-<compiler spec>.cmake``
   Give these files to your regular user account
   and log back in to that account. Copy these new host-config files to the
   ``host-configs`` subdirectory and commit them to your branch. Make sure all
   file changes from all previous steps are also committed and pushed upstream.

#. **Build new Docker images.**
   We utilize Docker images that have pre-built TPLs in our Github CI checks.
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

#. Make sure all changes in your branch are committed and pushed, and create
   a pull request for a merge to develop. If everything went well, all checks
   on your Github PR should pass.
 