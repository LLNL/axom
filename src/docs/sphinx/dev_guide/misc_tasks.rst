.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _misctasks-label:

********************************
Miscellaneous Development Items
********************************

This section describes various development tasks that need to be 
performed that are not covered in earlier sections.


===================
Web Documentation
===================

Describe how to build and install web documentation...

Shared LC web content location ``axom/src/docs/sphinx/web``


==================================
Third-party Library Installation
==================================

Describe how to run the scripts to install third-party libraries for 
testing different versions locally on a branch and for installing new
libraries for the team to use...

Building and installing TPLs for all compilers on the current LC platform you are on::

   $ python ./scripts/llnl_scripts/build_tpls.py -d <output/path/>

Questions we need to answer include:

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    and experimentation?
  * Others?

.. note :: Pull in content from ``../web/build_system/thirdparty_deps.rst`` ...
           fill in gaps and make sure it it up-to-date...

Updating a TPL to a new version
-------------------------------

Follow these steps to update Axom to use a new released version of one
of the TPL libraries.

#. Create a branch to work on the task.


#. Update the ``packages.yaml`` files in ``scripts/spack/configs/``.
   Find the entry for the library you wish to update and change the
   version number. Do this for the ``packages.yaml`` file of each system type,
   including those in the `docker` subdirectory.


#. Update ``scripts/uberenv/packages/[package name]/package.py`` for the
   library you are changing. Usually this can be done by copying directly
   from ``var/spack/repos/builtin/packages/[package name]/package.py`` in
   the Spack distribution. Clone the repo at ``github.com/spack/spack`` to
   find this file and copy it. Verify that the new version you intend to
   use is included in this file in the list of ``version()`` entries.

#. Do a test installation in your local directories

   $ scripts/llnl/build_tpl.py -d ../your/local/install/path

   This will do a full installation of the TPLs, which you can check to verify
   that the correct new version is installed. Also, it produces new host-config
   files that you should use to build and test Axom with this installation.
   These new host-config files will be located at the base of your local
   clone of the repository. If any changes to Axom code are needed to work
   with the new TPL update, make these changes and test them.

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
