.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _misctasks-label:

********************************
Miscellaneous Development Items
********************************

This section describes various development tasks that need to be 
performed at times and which are not covered in other sections.

We need to check that answers to the following questions are contained in
this section and are clear:

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    and experimentation?
  * Others?


===================
Web Documentation
===================

Axom web-based documentation is hosted on our 
`Read the Docs project <https://readthedocs.org/projects/axom/>`_. 
Multiple versions are visible there, including the latest content on the 
develop branch (*latest*) and the main branch (*main*). The documentation 
that appears is automatically re-generated each time a change is pushed to 
a branch in the GitHub repository that is enabled to be displayed on the 
Read the Docs project. If you are modifying Axom documentation, you can enable 
the branch you are working on so that you can see what it looks like as you 
push changes to the branch. If your documentation changes are part of a GitHub
pull request, it is a good idea to enable the documentation for that branch
and put a link to it in the pull request summary. This makes it easy for 
reviewers to check over your changes.

.. note :: When you no longer need the documentation of your branch to be
           visible on Read the Docs (e.g., your pull request is merged), 
           please disable that branch on Read the Docs.


========================================================================
Axom Third-party Library (TPL) Dependency Installation and Update
========================================================================

Axom dependencies are grouped into three categories: Git submodules,
system TPL dependencies, other TPL libraries. The following sections
describe how to install and update these dependencies for Axom.

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

TPLs
----

The following instructions describe how to install local copies of Axom
TPLs on Livermore Computing (LC) platforms. Typically, this process is 
followed when you want to build and test Axom against new versions of TPLs
to make sure everything works properly before deploying a new TPL set for
other Axom developers to use, to use in Axom CI testing, etc.

#. **Working on a local branch.** 
   Make a local clone of the Axom repo and create a branch to work on. Working
   on a separate branch is a good practice to avoid propagating mis-steps
   to other users and/or developers.

#. **Changing versions of system packages or other TPLs.**
   To change a version of a system package, which applies to an LC platforms 
   or a Docker container image we use for CI testing on GitHub, go into
   the directory ``axom/scripts/spack/configs``. There you will find a 
   sub-directory for each supported LC system type. Each sub-directory
   has a ``packages.yaml`` file which contains an entry for each system level
   package we rely on. Find the entry for the library you wish to update and 
   change the version number. Do this for each system you want to test/change,
   including configurations in the `docker` subdirectory.

   .. note :: Alongside each ``packages.yaml`` in each system package directory,
              there is a ``compilers.yaml`` file containing compiler and 
              version information for compilers we use for development and 
              testing. If you wish to test and build with a new compiler or 
              version on a system, modify the appropriate ``compilers.yaml`` 
              file.

   To change a version of a non-system TPL, go into the 
   ``axom/scripts/spack/packages`` directory. There you will find a 
   sub-directory for each TPL Axom uses. Modify the contents of the Spack
   package file ``package.py`` in each package sub-directory as needed. 

#. **Test installation.**
   Do a test installation of your changes in a local directory by running
   the following command in the top-level Axom directory::

   $ ./scripts/llnl_scripts/build_tpls.py -d local/install/path

   where ``local/install/path`` is a directory location where you want the 
   libraries to be installed.

   This will do a full installation of the TPLs, which you can use to check 
   that the the correct new version is installed and works. 

   Running the command above also produces new host-config files (i.e., 
   CMake cache files) that you use to build and test Axom with the installation.
   These new host-config files will be located in the top-level Axom directory
   of your local clone of the repo. If any changes to Axom code are needed 
   to work with the TPL update(s), make the changes and test them.

#. When you are confident that everything is correct, log in as user
   ``atk`` to each of the machines named in Axom's standard host-configs and run

   $ scripts/llnl/build_tpl.py

   This will do all of the standard installations in the shared directories
   used by Axom developers. When completed, they will produce new host-config
   files for each configuration. Give these files to your regular user account
   and log back in to that account. Copy these new host-config files to the
   ``host-configs`` subdirectory and commit them to your branch. Make sure all
   file changes from all previous steps are also committed and pushed upstream.

#. Next, build the docker images for continuous integration using GitHub
   actions. From Axom's GitHub page, click on "Actions" and then on "Docker
   TPL build" in the "Workflows" menu. Find the "Run Workflow" drop-down
   menu, select your branch, and click on the "Run workflow" button. This
   will launch the build of the docker images.

#. When the docker image build completes, click on your build and find the
   "Artifacts" listed at the bottom of the page. These contain host-configs
   for building Axom on the docker images. Download them and copy them to
   Axom's ``host-configs/docker`` subdirectory.

#. To complete the setup of the new docker images, the Compiler_ImageName
   entries in ``azure-pipelines.yaml`` at the top-level directory must be updated
   with the timestamped names of the new images. The new names can be found in
   the log files from the successful GitHub action. On the left of the page for
   the successful action is a "Jobs" menu. Click on each job and then find
   the "Build and push" section of the log. Within the first few lines of the
   section there should be an entry of the form
   ``"tags: axom/tpls:clang-10_12-18-20_00h-10m``. Copy the name beginning with
   ``axom/tpls`` to the appropriate locations in azure-pipelines.yaml. Repeat
   this with the names from each compiler job used in the GitHub action.

#. Make sure all changes in your branch are committed and pushed, and create
   a pull request for a merge to develop.
 

===================
Code Health Tools
===================

This section describes how to run code health tools we use.


Code Coverage
---------------

Setting up and running code coverage analysis...


Static Analysis
---------------

Setting up and running static analysis tools....


Memory Checking
----------------

Setting up and running memory checking tools....
