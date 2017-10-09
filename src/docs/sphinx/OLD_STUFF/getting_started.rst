.. ##
.. ## Copyright (c) 2017, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## LLNL-CODE-xxxxxx
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

================================
Getting Started with Axom
================================

This page gives details about getting started using Axom on LLNL LC platforms where third-party dependencies are already installed.
Note: Until we have a more formal release process in place, please follow the instructions below to use Axom.

---------------------
Checking out the code  
---------------------

The Axom source code can be obtained from the axom git/bitbucket project located at
  `<https://lc.llnl.gov/bitbucket/projects/ATK>`_

If you've never used git/bitbucket, you will need to create an `SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and  `add that SSH Key <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_ to your `CZ Bitbucket Profile <https://lc.llnl.gov/bitbucket/account>`_.

1.  First, clone the Axom repo in some local directory you own.::

     git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git

     cd axom

2.  Run the "setup-for-development" script, which will ensure that Git is configured correctly and client-side hooks are installed:::

            ./scripts/setup-for-development

3.  Then checkout the branch you wish to use. For example, the *"develop"* branch is currently the most robust and tested branch. Axom development work is pushed to it regularly.::

            git checkout develop

    *   If you are using a development branch of Axom and wish to keep up with the latest developments, you will need to pull them into your local copy.  This can be done as follows: ::

            git pull

    *  If you wish to follow and use a different Axom branch (i.e., a topic branch), you can check it out with the following commands: ::

         # update the list of remote branches
           git fetch

         # list the available branches

           git branch -a 

           git checkout {branch name, excluding the remotes/origin/}

   Here is a concrete example: ::

       git branch -a | grep unify 

           remotes/origin/feature/zagaris2/restructure-and-unify

       git checkout feature/zagaris2/restructure-and-unify

       Branch feature/zagaris2/restructure-and-unify set up to track remote branch feature/zagaris2/restructure-and-unify from origin.

       Switched to a new branch 'feature/zagaris2/restructure-and-unify'

If you wish to do development on Axom and push your work back to the Axom git/bitbucket repo, please see the instructions for `Code Development in the Toolkit <https://lc.llnl.gov/confluence/display/ASCT/Code+Development+in+the+Toolkit>`_.

-------- 
Git Tips
--------
1) Get your prompt to show which branch you are currently on, whether there are untracked (i.e., new) files, modified files, staged, unstaged or stashed changes. See `git-prompt.sh <https://github.com/git/git/blob/master/contrib/completion/git-prompt.sh>`_ 
2) Enable tab-autocompletion for git commands. See `git-completion <https://github.com/git/git/tree/master/contrib/completion>`_.

------------------------------------
Configuring and Building Axom
------------------------------------
The Axom build system requires that you use a version of CMake greater than 3.1.  You can download cmake binaries or source from:
 * `<http://www.cmake.org/download/>`_

Once you have a working version of CMake, here are two ways to configure and build axom:

------------------------------------------------------
Configuring using the 'configure' python helper script
------------------------------------------------------
If you have python on your platform, you can use the script "scripts/config-build.py" to setup a build. This script is designed to locate the correct cmake cache file to use from the 'host-config' directory, based on your platform and options you have provided.  Scripts options can be displayed by running 'config-build.py --help'.

**Example:** ::
 
 1.  From the root of the axom repo, with cmake in your path run:
       ./scripts/config-build.py

     # If no options are provided, will configure and create a build directory with the default compiler for this platform

 2.  cd into the build directory; e.g., 
       cd build-chaos-gnu-debug

 3.  Build axom's libs and unit tests
         make 

 4.  Run the axom unit tests
         make test

--------------------------------
Configuring manually using CMake
--------------------------------
You can also execute CMake directly to configure a build. ::

 1. From the root of the axom repo:
    # our CMake setup disallows in-source builds, so you need to create a build directory

         mkdir mybuild
         cd mybuild

 2. Configure an out-of-source build in a directory named "build-debug"
         cmake  -DCMAKE_BUILD_TYPE=Debug ../src

     or, to include a default cache file you can use the -C option

         cmake  -DCMAKE_BUILD_TYPE=Debug -C ../host-configs/other/Darwin.cmake ../src

         cmake  -DCMAKE_BUILD_TYPE=Debug -C ../host-configs/$SYS_TYPE.cmake ../src

 3. Build axom's libs and unit tests
        make 

 4. Run the axom unit tests
        make test

You can use ccmake or cmake-gui to modify specific build options.

* cmake:     `<https://cmake.org/cmake/help/v3.0/manual/ccmake.1.html>`_
* cmake-gui: `<https://cmake.org/cmake/help/v3.0/manual/cmake-gui.1.html>`_



