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

.. _release-label:

*******************************************
Axom Release Process
*******************************************

The Axom team uses its collective professional judgement to decide when
the code is ready for a release. Most often, a release is done when some
set of functionaly is sufficiently complete to make available to users or
when some other development goal is reached. This section describes how 
Axom releases are done. The process is fairly informal; however, it should 
be adequate to ensure that the software is in a reasonably robust and stable 
state when a release is done. In particular, we follow this process to avoid
simple oversights and issues that we do not want to pass on to users.

In the :ref:`gitflow-label` section, we noted that the **master branch
records the official release history of the project**. Specifically,
whenever, the master branch is changed, it is tagged with a new
version number. We use a git 'lightweight tag' for this purpose. Such
a tag is essentially a pointer to a specific commit on the branch.

Here are the steps to generate a new tagged version of the master branch
for release:

  #. Create a **release branch** off of the develop branch. Rather than
     finalize the release directly on the develop branch, a release branch
     is preferred so that other work may continue on the develop branch
     without interruption.

  #. Update the version number information in the
     ``axom/src/cmake/AxomVersion.cmake`` file.

  #. Create a pull request to merge the release branch into master so that
     it can be reviewed.

     * Ensure the code builds on all platforms and all tests pass.

     * Make sure any discovered issues are resolved.

     * Make sure all documentation (user guides, release notes, etc.) is
       updated and reviewed.

  #. When the release branch is ready and approved, merge it into master.

  #. Also merge master into develop if the release branch changed during
     this process.

  #. Tag the master branch with a new version number. See :ref:`semver-label`
     for a description of how version numbers are chosen. To tag master::

       $ git checkout master
       $ git tag vMM.mm.pp

     Recall that 'MM' is the major version number, 'mm' is the minor version
     number, and 'pp' is the patch number. The command::

       $ git tag

     will list all the tags and you should see the new tag you created.
     The command::

       $ git show vMM.mm.pp

     will show information about the tagged commit.

  #. Push the new tag to bitbucket.  You can push a single tag
     to bitbucket by running the command ::

       $ git push vMM.mm.pp

     To push all local tags to bitbucket, use this command::

       $ git push --tags

