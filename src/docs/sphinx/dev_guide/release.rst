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

The Axom team decides as a group when the code is ready for a release. 
Typically, a release is done when we want to make changes available to users;
e.g., when some new functionality is sufficiently complete or we want users to
try something out and give us feedback early in the development process. A
release may also be done when some other development goal is reached. This 
section describes how an Axom releases is done. The process is fairly 
informal. However, we want to ensure that the software is in a reasonably 
robust and stable state when a release is done. We follow this process to 
avoid simple oversights and issues that we do not want to pass on to users.

In the :ref:`gitflow-label` section, we noted that the **master branch
records the official release history of the project**. Specifically,
whenever, the master branch is changed, it is tagged with a new
version number. We use a git 'lightweight tag' for this purpose. Such
a tag is essentially a pointer to a specific commit on the master branch.

We finalize preparations for a release on a release branch so that other 
work may continue on the develop branch without interruption. 

.. note:: No significant code development is performed on a release branch.
          In addition to preparing release notes and other documentation, the
          only code changes that should be done are bug fixes identified 
          during release preparations

Here are the steps to follow for an Axom release.

  #  **Create a release branch** off of the develop branch to finalize 
     the release. The name of a release branch must contain the associated 
     release version number. Typically, we use a name like rc0.5.0 (i.e., 
     version 0.5.0 release candidate). See :ref:`semver-label` for a 
     description of how version numbers are chosen. 

  #. **Create a pull request to merge the release branch into master** so that
     release changes can be reviewed. Such changes include:

     * Update the version information (major, minor, and patch version numbers)
       at the top of the ``axom/src/cmake/AxomVersion.cmake`` file.

     * Test the code by running it through all continuous integration tests 
       and builds. This will ensure that all build configurations are working 
       properly and all tests pass. 

     * Fix any issues discovered during final release testing if code changes
       are reasonably small and re-run appropriate tests to ensure issues are
       resolved. If a major bug is discovered, and it requires significant
       code modifications to fix, do not fix it on the release branch. Add
       a JIRA issue for it and note the issue in the ``known bugs`` section of
       the release notes.

     * Make sure all documentation (source code, user guides, etc.) is
       updated and reviewed. This should not be a substantial undertaking as
       most of this should have been done during the regular development cycle.

     * Proofread the release notes for completeness and clarity and address
       any shortcomings. Again, this should not take much time as release notes
       should be updated during the regular development cycle. See 
       :ref:`release-notes-label` for information about release notes.

  #. **Merge the release branch into the master branch** when the release 
     branch is ready and approved.

  #. **Tag the master branch with a new version number.** To tag master::

       $ git checkout master
       $ git tag vMM.mm.pp

     Recall that 'MM' is the major version number, 'mm' is the minor version
     number, and 'pp' is the patch number. The command::

       $ git tag

     will list all the tags and you should see the new tag you created.
     The command::

       $ git show vMM.mm.pp

     will show information about the tagged commit.

  #. **Push the new tag to bitbucket.** You can push a single tag
     to bitbucket by running the command ::

       $ git push vMM.mm.pp

  #. **Create a pull request to merge master into develop.** When approved,
     merge it.


.. _release-notes-label:

*******************************************
Release Notes
*******************************************

Axom release notes are maintained in a single file ``axom/RELEASE-NOTES.md``.
The release notes for the latest version are at the top of the file. 
Notes for previous releases appear after that in descending version number 
order.

For each version, the release notes must contain the following information:

 * Axom version number and date of release

 * One or two sentence overview of release, including any major changes.
 
 * Release note items should be broken out into the following sections: 

    * Added: Descriptions of new features
    * Removed: Notable removed functionality
    * Deprecated: Deprecated features that will be removed in a future release
    * Changed: Enhancements or other changes to existing functionality
    * Fixed: Major bug fixes
    * Known bugs: Existing issues that are important for users to know about

.. note:: Release notes for each Axom version should explain what changed in 
          that version of the software -- and nothing else!! 

Release notes are an important way to communicate software changes to users
(functionality enhancements, new features, bug fixes, etc.). Arguably, they 
are the simplest and easiest way to do so. Each change listed in the release 
notes should contain a clear, concise statement of the change. Items should 
be ordered based on the impact to users (higher impact - first, lower impact 
last). 

.. note:: When writing release notes, think about what users need to know and 
          what is of value to them. 

Release notes should summarize new developments and provide enough detail
for users to get a clear sense of what's new. They should be brief -- don't
make them overly verbose or detailed. Provide enough description for users
to understand a change, but no more than necessary. In other words, release
notes summarize major closed issues in a human-readable narrative. Direct 
users to other documentation (user guides, software documentation, example 
codes) for details and additional information.

.. note:: Release notes should be updated as work is completed and reviewed
          along with other documentation in a pull request. This is much 
          easier than attempting to compile release notes before a release 
          by looking at commit logs, etc. **Preparing release notes as part 
          of the release process should take no more than one hour.

Lastly, release notes provide an easy-to-find retrospective record of 
progress for users and other stakeholders. They are useful for developers
and for project reporting and reviews.


