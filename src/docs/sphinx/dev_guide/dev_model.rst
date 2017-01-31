.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

****************************************
CS Toolkit Development Process Overview
****************************************

The CS Toolkit team uses an agile, sprint-based development process where
we target a set of development tasks for each release. We also employ a 
variety of tools to manage our workflow and software development. This guide
contains descriptions of our processes and information about how we use tools 
and about getting started with them.

The CS Toolkit Git repository lives in a
`Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_
on the Livermore Computing CZ Atlassian Space. The repository is the central
interaction hub for CS Toolkit code development.

=======================
First Things First...
=======================

The CS Toolkit team uses various means of communication and interaction.

The main interaction hub for developers is the Atlassian tool suite on the
Livermore Computing Collaboration Zone (CZ). These tools can be accessed through
the `MyLC Portal <https://lc.llnl.gov>`_.  Access to our Atlassian project 
spaces requires membership in one or both of the LC groups 'toolkit' and 
'toolkitd'. The first is for users and the second is for developers -- 
team members are in both. If you're not in these groups, send email to 
'asctoolkit-dev@llnl.gov' and request to be added.

We maintain two email lists: 'asctoolkit-dev@llnl.gov' 
(mainly for developer communication) and 'asctoolkit-users@llnl.gov' 
(for developer-user interaction). Go to the 
`LLNL E-Mail List Manager Site <https://lists.llnl.gov>`_ to subscribe to
these lists.

We also have a chat room on LLNL's Cisco Jabber instance called 
'CS Toolkit Dev'. It's open to anyone. You just have to log on to Jabber and 
join.


.. _gitflow-label:

======================================================
Gitflow Branching Model
======================================================

The CS Toolkit team follows the 'Gitflow' branch development model, which is
summarized in this section. See the `Atlassian Gitflow Description <https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow>`_ for more details.

Gitflow is a branching model centered around software releases. It is a simple 
workflow that makes clear which branches correspond to which phases of 
development and those phases are represented explicitly in the structure of 
the repository. As in other branching models, developers develop code locally 
and push their work to a central repository. 

----------------------------
Master and Develop Branches
----------------------------

The master and develop branches are the two main branches used in Gitflow.
They always exist and the distinction between them is central to the Gitflow
model. Other branches are temporary and used to perform specific development 
tasks.

The master branch records the official release history of the project. 
Each time the master branch is changed, it is tagged with a new version number.
For a description of our versioning scheme, see :ref:`versioning-label`.

The develop branch is used to integrate and test new features and most 
bug fixes before they are merged into master. 

.. important:: **Development never occurs directly on the master or develop
               branches.**  

---------------
Topic Branches
---------------

Topic branches are created off of other branches (usually develop)
and are used to develop new features and resolve issues before they 
propagate to master. Topic branches are temporary, living only as long as they
are needed to complete a development task.

Each new feature, or other well-defined portion of work, is developed on its 
own topic branch, with changes being pushed to the central repository regularly
for backup. We typically include a label, such as  "feature" or "bugfix", in 
the topic branch name to make it clear what type of work is being done on the 
branch. See :ref:`topicdev-label` for a description of common Git mechanics 
when doing topic branch development. 

When a feature is complete, a pull request is submitted for review by other 
team members. When all issues arising in a review have been addressed and 
reviewers have approved the pull request, the feature branch is merged into 
develop. See :ref:`review-label` for more information about code reviews and 
pull request approval.

.. important:: **Feature branches never interact directly with the master 
               branch.**

-----------------
Release Branches
-----------------

Release branches are another important temporary branch type in Gitflow:
When the team has decided that enough features, bug fixes, etc. have been 
merged into develop (for example, all items identified for a release have 
been completed), a release branch is created off of develop to finalize the 
release. Creating a release branch starts the next release cycle on develop. 
At that point, new work can start on feature branches for the next release. 
Only changes required to complete the release are added to a release branch. 
When a release branch is ready, it is merged into master and master is tagged 
with a new version number. Finally, master is merged back into develop since 
it may have changed since the release was initiated.

.. important:: **No new features are added to a release branch. Only bug fixes, 
               documentation, and other release-oriented changes go into a 
               release branch.**

----------------
Hotfix Branches
----------------

The last important temporary branch type in Gitflow is a hotfix branch.
Sometimes, there is a need to resolve an issue in a released version on the 
master branch. When the fix is complete, it is reviewed using a pull request 
and then merged into both master and develop when approved. At this point, 
master is tagged with a new version number. A dedicated line of development 
for a bug fix, using a hotfix branch, allows the team to quickly address 
issues without disrupting other parts of the workflow. 

.. important:: Hotfix branches are the only branches created off of master.

-----------------------
Gitflow Illustrated
-----------------------

The figure below shows how branches interact in Gitflow.

.. figure:: gitflow-workflow.png

   This figure shows typical interactions between branches in the Gitflow 
   workflow. Here, master was merged into develop after tagging version v0.1. 
   A fix was needed and so a hotfix branch was created. When the fix was 
   completed, it was merged into master and develop. Master was tagged 
   with version v0.2. Also, work was performed on two feature branches. 
   When one feature branch was done, it was merged into develop. Then, a 
   release branch was created and it was merged into master when the release 
   was finalized. Finally, master was tagged with version v1.0.


.. _versioning-label:

======================================================
Semantic Versioning
======================================================

The CS Toolkit team follows the *semantic* versioning scheme for assigning
release numbers. Semantic versioning is a methodology for assigning version 
numbers to software releases in a way that conveys specific meaning about 
the code and modifications from version to version. 
See `Semantic Versioning <http://semver.org>`_ for a more detailed description.

----------------------------
Version Numbers and Meaning
----------------------------

Semantic versioning is based on a three part version number `MM.mm.pp`:

  * `MM` is the *major* version number. It is incremented when an incompatible 
    API change is made. That is, the API changes in a way that may break code
    using an earlier release of the software with a smaller major version 
    number. Following Gitflow (above), the major version number may be changed
    when the develop branch is merged into the master branch.
  * `mm` is the *minor* version number. It changes when functionality is
    added that is backward-compatible. The API may grow to support new 
    functionality. However, the software will function the same as any
    earlier release of the software with a smaller minor version number
    when used through the intersection of two APIs. Following Gitflow (above), 
    the minor version number is always changed when the develop branch is 
    merged into the master branch, except possibly when the major version 
    is changed.
  * `pp` is the *patch* version number. It changes when a bug fix is made that
    is backward compatible. That is, such a bug fix is an internal 
    implementation change that fixes incorrect behavior. Following Gitflow 
    (above), the patch version number is always changed when a hotfix branch
    is merged into master, or when develop is merged into master and the 
    changes only contain bug fixes.

-------------------------
Version Number Changes
-------------------------

A key consideration in meaning for these three version numbers is that
the software has a public API. Changes to the API or code functionality
are communicated by the way the version number is incremented. Some important
conventions followed when using semantic versioning are:

  * Once a version of the software is released, the contents of the release 
    *must not* change. If the software is modified, it *must* be released
    as a new version.
  * A major version number of zero (i.e., `0.mm.pp`) is considered initial 
    development where anything may change. The API is not considered stable.
  * Version `1.0.0` defines the first stable public API. Version number 
    increments beyond this point depend on how the public API changes.
  * When the software is changed so that any API functionality becomes 
    deprecated, the minor version number *must* be incremented.
  * A pre-release version may be denoted by appending a hyphen and a series
    of dot-separated identifiers after the patch version. For example,
    `1.0.1-alpha`, `1.0.1-alpha.1`, `1.0.2-0.2.5`.
  * Versions are compared using precedence that is calculated by separating
    major, minor, patch, and pre-release identifiers in that order. Major, 
    minor, and patch numbers are compared numerically from left to right. For 
    example, 1.0.0 < 2.0.0 < 2.1.0 < 2.1.1. When major, minor, and patch
    numbers are equal, a pre-release version has lower precedence. For 
    example, 1.0.0-alpha < 1.0.0.

By following these conventions, it is fairly easy to communicate intent of
version changes to users and it should be straightforward for users
to manage dependencies on the CS Toolkit.


.. _review-label:

======================================================
Code Reviews and Acceptance
======================================================

Before any code is merged into the develop or master branches, it
must be tested, reviewed, and accepted. Submitting a pull request on
the Toolkit Bitbucket project to merge a branch into develop or master 
initiates the test and review processes. All builds and tests must pass 
for a pull request to be approved. Also, it is expected that unit tests 
are constructed to exercise any new functionality that is introduced. This 
will be assessed by reviewers of each pull request. See :ref:`testing-label` 
for more information about testing.

Code changes in a pull request must be accepted by at least one member
of the Toolkit development team other than the originator of the pull
request. It is recommended to have several team members review pull 
requests, especially when changes affect APIs. Pull request reviewers can be 
selected on Bitbucket when the pull request is created. Changes reviewed by 
the team are accepted, rejected, or commented on for improvement; e.g., 
issues to be addressed, suggested changes, etc. Pull requests can be updated
with additional changes as needed. When a pull request is approved, it can 
be merged. If the merged branch is no longer needed for development, it 
should be deleted.

In addition to successful compilation and passing tests, changes to the 
develop and master branches should be scrutinized in other ways and using 
other tools. In particular :

* The code should compile cleanly at the highest warning level with the 
  main compilers supported by the project. All warnings **must** be 
  understood and eliminated if possible. Reducing a compiler warning 
  level to eliminate warning messages **is not** acceptable.

  Compiler warnings, while seemingly innocuous at times, often indicate
  problems that do not appear until later or until specific run-time
  conditions are encountered.

* Static analysis tools **should** be applied to the code using tools such
  as `cppcheck`, etc. to identify potential implementation issues.

* Runtime memory checking, using a  tool such as Valgrind, **should** be 
  performed to verify that there are no leaks or other memory issues. 

.. note :: Bamboo setup needs to be completed to automate builds and 
           testing for pull requests. We have not yet established policies 
           or included use of these tools in our Bamboo test plans. Ideally, 
           we would like to automate them as part of our CI and pull request 
           approval processes.


---------------------
Pull Request Summary
---------------------

  #. When code is ready to be considered for acceptance, create a pull request
     on the CS Toolkit Bitbucket project. Identify the appropriate reviewers 
     and add them to the pull request.

  #. Code must build successfully and all relevant tests must pass, including
     new tests required for new functionality.

  #. All issues (build failures, test failures, reviewer requests) must be 
     addressed before a pull request is accepted.

  #. Pull requests must be approved by at least one member of development 
     team other than the pull request originator.

  #. When a pull request is approved it may be merged. If the merged branch is
     no longer needed, it should be deleted. This can be done when merging
     with Bitbucket. 


---------------------------
Code Review Checklist
---------------------------

Beyond build and test correctness, we also want to ensure that code follows
common conventions before acceptance. The following list summarizes concerns 
we want to identify during pull request reviews and resolve before a pull 
request is approved for merging. Please see the CS Toolkit coding guidelines
**(insert link here)**
document for details. 

 #. A new file or directory must be located in in the proper location; e.g.,
    in the same directory with existing files supporting related functionality.
 #. File contents must be organized clearly and structure must be consistent 
    with conventions. 
 #. Namespace and other scoping conventions must be followed. 
 #. Names (files, types, methods, variables, etc.) must be clear, easily
    understood by others, and consistent with usage in other parts of the code.
    Terminology must be constrained; i.e., don't introduce a new term for 
    something that already exists and don't use the same term for different 
    concepts.
 #. Documentation must be clear and follow conventions. Minimal, but adequate, 
    documentation is preferred.
 #. Implementations must be correct, robust, portable, and understandable to
    other developers.
 #. Adequate tests (unit and performance) tests must be added for new 
    functionality.


