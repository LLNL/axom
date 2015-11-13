Getting Started with the Toolkit
---------------------------------

This page gives details about getting started using the Toolkit on LLNL LC platforms where third-party dependencies are already installed.
Note: Until we have a more formal release process in place, please follow the instructions below to use the Toolkit.
Checking out the code  
The Toolkit source code can be obtained from the Toolkit git/Stash project located at https://lc.llnl.gov/stash/projects/ATK. If you've never used git/Stash, you will need to create an SSH key and add that SSH Key to your CZ Stash Profile.
First, clone the Toolkit repo in some local directory you own.
git clone ssh://git@cz-stash.llnl.gov:7999/atk/asctoolkit.git
cd asctoolkit
Then checkout the branch you wish to use. For example, the "develop" branch is currently the most robust and tested branch. Toolkit development work is pushed to it regularly.
git checkout develop
If you are using a development branch of the Toolkit and wish to keep up with the latest developments, you will need to pull them into your local copy.  This can be done as follows:
git pull
If you wish to follow and use a different Toolkit branch (i.e., a topic branch), you can check it out with the following commands:
# update the list of remote branches
git fetch
# list the available branches
git branch -a 
git checkout {branch name, excluding the remotes/origin/}
Here is a concrete example:
git branch -a | grep unify 
   remotes/origin/feature/zagaris2/restructure-and-unify
git checkout feature/zagaris2/restructure-and-unify
Branch feature/zagaris2/restructure-and-unify set up to track remote branch feature/zagaris2/restructure-and-unify from origin.
 Switched to a new branch 'feature/zagaris2/restructure-and-unify'

If you wish to do development on the Toolkit and push your work back to the Toolkit git/Stash repo, please see the instructions for Code Development in the Toolkit.
Git Tips
1) Get your prompt to show which branch you are currently on, whether there are untracked (i.e., new) files, modified files, staged, unstaged or stashed changes. See git-prompt.sh 
2) Enable tab-autocompletion for git commands. See git-completion.
Configuring and Building the Toolkit
The Toolkit build system requires that you use a version of CMake greater than 3.1.  You can download cmake binaries or source from: http://www.cmake.org/download/
Once you have a working version of CMake, here are two ways to configure and build the toolkit:
Configuring using the 'configure' python helper script
If you have python on your platform, you can use the script "scripts/config-build.py" to setup a build. This script is designed to locate the correct cmake cache file to use from the 'host-config' directory, based on your platform and options you have provided.  Scripts options can be displayed by running 'config-build.py --help'.
Example: 
# From the root of the asctoolkit repo, with cmake in your path run:
./scripts/config-build.py
# If no options are provided, will configure and create a build directory with the default compiler for this platform
# cd into the build directory; e.g., 
cd build-chaos-gnu-debug
# build the toolkit libs and unit tests
make 
# run the toolkit unit tests
make test

Configuring manually using CMake
You can also execute CMake directly to configure a build.
# from the root of the asctoolkit repo
# our CMake setup disallows in-source builds, so you need to create a build directory
mkdir mybuild
cd mybuild
# configure an out-of-source build in a directory named "build-debug"
cmake  -DCMAKE_BUILD_TYPE=Debug ../src
#
# or, to include a default cache file you can use the -C option
# cmake  -DCMAKE_BUILD_TYPE=Debug -C ../host-configs/other/Darwin.cmake ../src
# cmake  -DCMAKE_BUILD_TYPE=Debug -C ../host-configs/$SYS_TYPE.cmake ../src
#
# build the toolkit libs and unit tests
make 
# run the toolkit unit tests
make test

You can use ccmake or cmake-gui to modify specific build options.
TODO: George will look for links to docs from kitware about ccmake and cmake-gui.


