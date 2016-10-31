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

.. _tooleco-label:

======================================================
Software Development Tool Ecosystem
======================================================

We use a variety of tools for software development. Our tool philosophy has
three main tenets:

  * Use robust, commonly-used tools - don't invent something if we don't need to
  * Use the tools in ways that are easy for non-experts
  * Strive for automation and reproducibility

The main tools we use are listed in this section. Details about how we use 
them and helpful information about getting started are provided 
in the sections that follow.

* We use our `Confluence project space <https://lc.llnl.gov/confluence/display/ASCT/ASC+Simulation+CS+Toolkit+Home>`_ for team discussion, planning, maintaining meeting notes, etc.
* Our Git repository houses the Toolkit source code, build configurations, scripts, test suites, documentation, etc. The repository lives in our `Bitbucket project <https://https://lc.llnl.gov/bitbucket/projects/ATK>`_
* We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for issue tracking.
* We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for continuous integration and automated testing.

---------------------
Build, Link, Triumph
---------------------

Our build system, called *BLT*, is maintained in its own repo in our 
Bitbucket project. **Add link to BLT documentation when it is available** 
BLT provides a "common sense" setup based on CMake. It manages our build 
environment (compilers, programming models - OpenMP, MPI, CUDA, etc., and 
third-party library locations) as well as our software development tool 
integration via *make targets*. BLT has built-in support for the following
tools, all of which we use for Toolkit development:

* Documentation - *Doxygen* (source code docs) and *Sphinx* (user docs)
* Unit testing - *CTest* (test orchestration), *Google Test* (C/C++ unit tests), *Fruit* (Fortran unit tests)
* Code Health - *Uncrustify* (code style), *gcov* and *lcov* (code coverage), and *Valgrind* (memory checking)
* Benchmarking - *Google Benchmark*

We use `Spack <https:://github.com/LLNL/spack>`_ to manage and build the 
third-party libraries on which the Toolkit depends. We also maintain some
custom python scripts (in our 'scripts directory') to bootstrap Spack.

More information on building the code and third-party libraries can be found
in :ref:`configbuild-label`.

.. note :: BLT is supported as a standalone product and used by other 
           software projects.


--------------------------------------
Git/Bitbucket
--------------------------------------

This section provides some information about getting started with Git and 
Bitbucket and describes operations related to topic branch development 
on the CS Toolkit project. Our Git repository lives in our 
`Bitbucket project <https://https://lc.llnl.gov/bitbucket/projects/ATK>`_.

If you are new to the Git or want to brush up on its features, the 
`Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of good 
information, and `Learn Git Branching <http://learngitbranching.js.org/>`_ 
is nice for visual, hands-on learners. Also, the book 
`Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an
excellent comprehensive guide to Git. 

To make Git easier to work with, folks have written some useful scripts. For 
example, see `Git scripts <https://github.com/git/git/tree/master/contrib/completion>`_ for scripts that enable tab-autocompletion for Git commands and set 
up you prompt to show which branch you are on, etc.

SSH keys
^^^^^^^^^^^^^^^

If you have not used Bitbucket before, you will need to 
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

.. _repoclone-label:

Getting a local working copy of the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before starting work on the code, you must clone the repo into your working
space. This is done by typing::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/asctoolkit.git

Important notes:

  * You don't need to remember the URL for the Toolkit repo above. It can be 
    found by going to the CS Toolkit repo on our Bitbucket project and 
    clicking on the 'Clone' Action button that appears when you hover your 
    mouse cursor over the ellipses on the top left of the web page.
  * The '--recursive' argument above is needed to pull the BLT build system into
    your local copy of the repo. It is a Git sub-module of the Toolkit.

After cloning, enter the top-level Toolkit directory and run the development
setup script we provide to ensure that your Git environment is configured 
properly and client-side hooks we use are installed; i.e.,::

  $ cd asctoolkit
  $ ./scripts/setup-for-development.sh

More about the Git hooks later.

.. _topicdev-label:

Performing topic branch development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is worth re-emphasizing a fundamental fundamental principle of the Gitflow 
development model that we described in :ref:`gitflow-label`.

.. important:: **We never work directly on the develop or master branches. 
               We use topic branches instead.**

When we refer to a *topic branch*, it could be a *feature branch*, 
a *bugfix branch*, etc. The basic workflow for performing development 
on a topic branch is:

  #. Create a topic branch off the develop branch and push the new branch
     to Bitbucket.
  #. Make changes and commit them to your branch in your local copy of the repo.
     Remember to push changes to the main repo on Bitbucket regularly for 
     backup.
  #. If you are working on your topic branch for a while, it is a good idea
     to keep your topic branch current with develop by merging develop into
     your topic branch regularly.
  #. When your work is complete, create a pull request so others on the team 
     can review your work. See :ref:`review-label`.

Here are some details about each of these steps.

**Step 1.** A topic branch name should include your user id and a brief 
description indicating the purpose of the branch. We typically label such 
branches using "feature", "bugfix", etc. to make it clear what type of work 
is being performed on a branch. For example,::

  $ git checkout -b feature/<userid>/some_cool_new_feature
  $ git push -u

You can also attach a JIRA issue number to the branch name. Then, Bitbucket 
will associate the issue with the commit when you merge your branch to the 
develop branch. For example,::

  $ git checkout -b bugfix/<userid>/jira-atk-<issue #>
  $ git push -u

In each of these examples, the 'git push -u' command pushes the branch to 
Bitbucket and it will appear in the list of branches you and other developers 
can see there.

**Step 2.** After the topic branch is created, and you've pushed it to 
Bitbucket, perform your development; i.e., edit files, add files, etc. 
Common commands you will use are::

  $ git add <file>
  $ git commit
  $ git push 

The 'add' command adds a file (or files) to be staged for commit. The 'commit'
command commits staged files to your local copy of the repository. The 'push' 
command pushes your commits to the topic branch in the main Git repo. You 
could also do::

  $ git push origin

This is equivalent to 'git push' if you specified the '-u' option when you
originally pushed your topic branch you created it.

Recall the Git environment setup script we recommended that you run after
cloning the repo in the Section :ref:`repoclone-label` above. One of the
Git pre-commit hooks that the script sets up applies formatting constraints 
on the commit message you provide when you execute the 'commit' command. The
constraints are recommended Git practices that help make it easier to use 
various tools with the Git version control system.

**Step 3.** If you will be working on your branch for a while, it is a good 
idea to merge from the develop branch to your topic branch at reasonable 
intervals to avoid getting too far out of sync. Otherwise, you may have 
many conflicts to resolve when you are ready to merge your topic branch
into the develop branch and the merge could be difficult. 

Before you begin the merge, make sure all outstanding changes to your topic 
branch are committed. Then, you need to make sure your local repo is 
up-to-date with the main develop branch by checking it out and pulling in 
the latest changes; i.e.,::

  $ git checkout develop
  $ git pull

Next, you need to go back to your topic branch, merge changes in from the 
develop branch, and check for conflicts::

  $ git checkout <your topic branch>
  $ git merge develop

The 'merge' command will tell you whether there are conflicts and which
files have them. Hopefully, you will not see any conflicts and you can 
continue working on your topic branch. If there are conflicts, you must
resolve them before you will be able to merge your topic branch to develop.
So, you may as well resolve them right away. You can resolve them by
editing the conflicting files and committing the changes. Merge conflicts
appear in a file surrounded by lines with special characters on them. For
example, if you open a conflicted file in an editor, you may see::

  <<<<<<< HEAD
  // lines of code, etc...
  =======
  // more lines of code, etc...
  >>>>>>> develop

The first section is the file contents in current branch head (your topic 
branch). The second section is the version in the develop branch. To resolve
the conflict, choose the correct version of contents you want and delete the
other lines. 

Alternatively, you can use a tool to help resolve your conflicts. There is 
the 'git mergetool' command and the "meld" tool, which is very powerful and 
intuitive). Diff tools like "tkdiff" are also helpful for resolving merge 
conflicts.

.. important:: **Git will not let you commit a file with merge conflicts.**


**Step 4.** When you are ready to merge your topic branch to the develop 
branch, you must initiate a pull request in Bitbucket. This is done by going 
into the Toolkit Bitbucket project, selecting your branch, and clicking the 
pull request button -- make sure you select the correct destination branch. 
The default destination branch in our project is set up to be the develop 
branch. So, in most cases, you won't have to do anything special. 

You must also select appropriate team members to review changes. Our Bitbucket 
project is set up to require at least one other developer to approve the pull 
request.

.. important:: **You cannot approve your own pull request.**

When your pull request is approved (see :ref:`review-label` for more 
information), you merge your changes to the develop branch by clicking the 
"merge" button in Bitbucket. If there are no merge conflicts, the merge will 
proceed and you are done. If there are conflict, Bitbucket will tell you
before you try to merge.

If there are conflicts, Bitbucket will not allow the merge to proceed. 
You must resolve the conflicts first. The preferred way to do this is to go 
into your branch and do the following::

  $ git fetch origin
  $ git merge origin

The 'fetch' command pulls changes from the remote branch into your local 
branch. Running the 'merge' command will show which files have conflicts 
as we described in the previous step. Fix the conflicts as described in 
the previous step. After all conflicts are resolved, run the 'commit' and 
'push' commands as usual. 

  $ git commit
  $ git push

Lastly, complete the merge in Bitbucket by clicking the merge button.

.. important:: **To keep things tidy, please delete your topic branch in 
               Bitbucket after it is merged if you no longer need it for 
               further development. Bitbucket also provides an option to 
               do this before doing the merge.**

Checking out an existing branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When working on multiple branches, or working on one with someone else, you
will need to checkout a specific branch. Any existing branch can be checked 
out from the Git repository and cloned from, etc. Here are some useful 
commands::

  $ git fetch
  $ git branch -a
  $ git checkout <branch name>

The 'fetch' command updates the list of remote branches and the 'branch'
command lists the available branches. The 'checkout' command checks out
the specified branch into your local working space. **Note that you do not
give the '-b' option when checking out an existing branch.** The option is
only used when creating a new branch.

Here is a concrete example::
  
  $ git branch -a | grep homer
    remotes/origin/feature/homer/pick-up-bart
  $ git checkout feature/homer/pick-up-bart
    Branch feature/homer/pick-up-bart set up to track remote branch feature/homer/pick-up-bart
    Switched to a new branch 'feature/homer/pick-up-bart'


--------------------------------------
JIRA Issue Workflow
--------------------------------------

We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for 
issue tracking. There you can create issues, edit them, comment on them,
check issue status, and search for issues in various ways, including setting 
up filters to customize your searches.

Issue states
^^^^^^^^^^^^^

We have customized our issue workflow to make it simple and easy to understand.
Specifically, each issue has three possible states:

  * **Open.** Every issues starts out as an open issue. An open issue can 
    be assigned to someone or unassigned. When an issue is assigned, this 
    means that the assignee owns the issue and is responsible for working 
    on it. An open issue that is unassigned has not been been discussed or 
    reviewed, or no decision to act on it has been made. Typically, an open 
    issue means that it is not being worked on.
  * **In Progress.** An issue in progress is one that is actively being 
    worked on.
  * **Closed.** When an issue is closed, work on it has been completed, or 
    a decision has been made that it will not be addressed.

An open issue can transition to either in progress (work has started on it)
or closed. An in progress issue can transition to either open (work on it
has stopped, but it is not finished) or closed. Finally, a closed issue
can be re-opened, which changes its state to open.

The figure below shows issue state transitions in our JIRA workflow.

.. figure:: jira-issue.png

   This figure shows allowed state transitions in our JIRA issue workflow.


Creating a new issue
^^^^^^^^^^^^^^^^^^^^^

To create a new issue, click the 'Create' button at the top of the CS Toolkit
JIRA project page and enter information in the issue fields. Filling in the
fileds properly greatly helps other team members search through project issues
to find what they are looking for. Note that issue fields marked with a red 
asterisk are required. The others are not required, but may be used to include 
helpful information. The main issues we use regularly are:

  * **Project.** The CS Toolkit will show up as the default. You shouldn't need
    to change this.
  * **Issue Type.** We use only three issue types: *Bug*, *New Feature*, and
    *Task*. A bug is something broken that needs to be fixed. A new feature
    is something to add that increases functionality, enhances an interface,
    etc. Task is a "catch-all" issue type for any other issue.
  * **Summary.** Provide a short descriptive summary. A good (and brief)
    summary makes it easy to scan a list of issues to find one you are 
    looking for.
  * **Priority.** Select an appropriate issue priority to impart its level 
    of importance or urgency. Clicking on the question mark to the right of
    the priority field provides a description of each option.
  * **Components.** Each issue is labeled with the Toolkit component it 
    applies to. Other "component" labels indicate build system issues, 
    documentation issues, etc. 
  * **Assignee.** Unless you are certain which team member should be assigned
    the issue, choose 'Unassigned'. This will indicate that the issue requires
    discussion and review before it is assigned. The default assignee is the
    owner of the component you chose earlier if you make no choice.
  * **Reporter.** Unless you explicitly enter someone in this field, you, as
    the issue creator, will be the reporter. This is the correct choice in
    almost all cases.
  * **Description.** The description field should be used to include important
    details about the issue that will help the developer who will work on it.

Other fields that appear may be used also if you think they will help
describe the issue. However, the team seldom uses fields apart from the list
above.

Starting and stopping work on an issue
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you begin work on an issue, you should open it, make sure it is 
assigned to you and click the 'Start Progress' button at the top of the issue.
This changes its status to *In progress*.

If there is still work to do on the issue, but you will stop working on it 
for a while, you can click the 'Stop Progress' button at the top of the
issue. This changes its status back to *Open*.

Closing an issue
^^^^^^^^^^^^^^^^^

When work is completed on an issue (which includes testing, adding
new documentation if needed, etc.), or the issue will not be addressed,
it should be closed. To close an issue, click the 'Close' button and select 
the appropriate issue resolution. There are two options: *Done* and *Won't Fix*.
'Done' means that the issue is resolved. 'Won't Fix' means that the issue will 
not be addressed for some reason.

When closing an issue, adding information to the 'Comment' field may be 
helpful. For example, when an issue is closed as 'Won't Fix', it is helpful to
enter a brief explanation as to why this is so.

.. note :: BLT is supported as a standalone product and used by other 
           software projects.

Issue assignee
^^^^^^^^^^^^^^^

Note that an assigned issue can be assigned to someone else to work on it.
An assigned issue can also be set back to 'Unassigned' if it needs further 
discussion by the team.


--------------------------------------
Bamboo Continuous Integration
--------------------------------------

We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for 
continuous integration and automated testing. We maintain a collection of
test plans for performing automated and manual builds, tests, and other
code health monitoring tasks.

Automated plans
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note :: **Fill in this section with a description of these plans: what is
           built, tested, other tasks performed, when they are run, etc.**

Manually running a plan on a branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note :: **Fill in this section with a description of what can be run 
           manually, and how to do it.**

Restricted Zone (RZ) Bamboo Project
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note :: **Fill in this section with a description of this when it is
           set up.**

