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

The CS Toolkit team uses a variety of tools for software development. Our
guiding principles with respect to such tools are:

  * Adopt robust, commonly-used tools - don't invent something if we don't need to
  * Make the tools easy to use for non-experts
  * Strive for automation and reporoducibility

The main tools we use are listed in this section. Details about how we use 
these tools and helpful information about getting started are provided 
in the sections that follow.

The main interaction hub for developers is the Atlassian tool suite on the 
Livermore Computing Collaboration Zone (CZ). Access to our Atlassian project
spaces requires membership in the LC groups 'toolkit' and 'toolkitd'. The 
following list contains links to each of these spaces: 

* Our Git repository houses the Toolkit source code, build configurations, scipts, test suites, documentation, etc. The repository lives in our `Bitbucket project <https://https://lc.llnl.gov/bitbucket/projects/ATK>`_
* We use our `Confluence project space <https://lc.llnl.gov/confluence/display/ASCT/ASC+Simulation+CS+Toolkit+Home>`_ for team discussion, planning, maintaining meeting notes, etc.
* We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for issue tracking.
* We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for continuous integration and automated testing.

Our build system, called *BLT*, is maintained in its own repo in our 
Bitbucket project. **Add link to BLT documentation when it is available** 
BLT provides a "common sense" setup based on CMake. It manages our build 
environment (compilers, programming models - OpenMP, MPI, CUDA, etc., and 
third-party library locations) as well as our software development tool 
integration via *make targets*. BLT has built-in support for:

* Documentation - *Doxygen* (source code docs) and *Sphinx* (user docs)
* Unit testing - *CTest* (test orchestration), *Google Test* (C/C++ unit tests), *Fruit* (Fortran unit tests)
* Code Health - *Uncrustify* (code style), *gcov* and *lcov* (code coverage), and *Valgrind* (memory checking)
* Benchmarking - *Google Benchmark*

We use `Spack <https:://github.com/LLNL/spack>`_ to manage and build the 
third-party libraries on which the Toolkit depends and a custom python 
script to bootstrap Spack.

.. note :: BLT is supported as a standalone product and used by other 
           software projects.


--------------------------------------
Git/Bitbucket Basics
--------------------------------------

This section provides some information about getting started with Git and 
Bitbucket and describes operations related to topic branch development 
on the CS Toolkit project. 

Git resources
^^^^^^^^^^^^^^^

If you are new to the Git or want to brush up on some of its features, the 
`Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of good 
information, and `Learn Git Branching <http://learngitbranching.js.org/>`_ 
is nice for visual, hands-on learners. Also, the book 
`Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an
excellent comprehensive guide to Git. To make Git easier to work with, folks 
have written some useful scripts. For example, see `Git scripts <https://github.com/git/git/tree/master/contrib/completion>`_ for scripts that enable 
tab-autocompletion for Git commands and set up you prompt to show which branch
you are on, etc.

SSH keys
^^^^^^^^^^^^^^^

If you've not used Bitbucket before, you will need to 
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.

.. _repoclone-label:

Getting a local copy of the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Before starting work on the code, you must clone the repo into your working
space. This is done by typing::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/asctoolkit.git

Notes:

  * The URL above can be found by going to the CS Toolkit repo on our 
    Bitbucket project and clicking on the 'Clone' Action button that appears 
    when you hover your mouse cursor over the ellipses on the top left of 
    the web page.
  * The '--recursive' argument is needed to pull the BLT build system into
    your local copy of the repo. It is a Git submodule of the Toolkit.

After cloning, enter the top-level Tookit directory and run the development
setup script we provide to ensure that your Git environment is configured 
properly and client-side hooks are installed; i.e.,::

  $ cd asctoolkit
  $ ./scripts/setup-for-development.sh

More about the Git hooks later.

Performing topic branch development
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is worth emphasizing a fundamental principle of the Gitflow development
model that we described in :ref:`gitflow-label`:

.. important:: **We never work directly on the develop or master branches. We use topic branches instead.**

When we refer to a *topic branch*, it could be a *feature branch*, 
a *bugfix branch*, etc. The basic workflow for performing development 
on a topic branch is:

  #. Create a topic branch off the develop branch and push the new branch
     to Bitbucket.
  #. Make changes and commit them to your branch in your local copy of the repo.
     Push changes on your topic branch to Bitbucket regularly for backup.
  #. If you are working on your topic branch for a while, it is a good idea
     to keep your topic branch current with develop
  #. When your work is complete, create a pull request so others on the team 
     can review your work. See :ref:`review-label`.

Here are some details about these steps.

**Step 1.** A topic branch should include your user id and a description 
indicating the purpose of the branch. We typically label such branches using 
"feature", "bugfix", etc. to make it clear what work is being done on the 
branch. For example,::

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
Bitbucket, do your development in it; i.e., edit files, add files, etc. 
Common commands you will use are::

  $ git add <file>
  $ git commit
  $ git push 

The 'add' command adds a file (or files) to be staged for commit. The 'commit'
command commits staged files to your branch on your local copy of the 
repository. The 'push' command pushes your commits to the topic branch in 
the main Git repo. You could also do::

  $ git push origin

This is equivalenet to 'git push' if you specified the '-u' option to the
push command after you created your topic branch.

Recall the Git environment setup script we suggested that you run after
cloning the repo in the Section :ref:`repoclone-label` above. One of the
Git pre-commit hooks that it sets up applies formatting constraints on the
commit message you provide when you execute the 'commit' command. The
constraints are recommended Git practices that help make it easier to use 
various tools with the Git version control system.

**Step 3.** If you will be working on your branch for a while, it is a good 
idea to merge from the develop branch to your topic branch at reasonable 
intervals to avoid getting too out of sync. Otherwise, you may have to 
resolve many conflicts when you are ready to merge your topic branch work 
onto the develop branch.

Before you begin, make sure all of outstanding changes to your topic branch 
are committed. Then, you need to make sure your local repo is up-to-date with 
the main develop branch by checking it out and pulling in the latest changes. 
This can be done as follows::

  $ git checkout develop
  $ git pull

Next, you need to go back to your topic branch, merge changes in from the 
develop branch, and check for conflicts::

  $ git checkout <your topic branch>
  $ git merge develop
  $ git status

The 'status' command will tell you whether there are merge conflicts and which
files have them. Hopefully, you will not see any conflicts and you can 
continue working on your topic branch. If there are conflicts, you must
resolve them before you will be able to merge your topic branch to develop.
So, you may as well resolve them right away. You can resolve them by
editing the conflicting files and committing the changes or you can use a 
tool to help resolve your conflicts. There is the 'git mergetool' command and
the "meld" tool (very powerful and intuitive). Diff tools like "tkdiff" are 
also helpful for resolving merge conflicts.

**Step 4.** When you are ready to merge your topic branch to the develop 
branch, you must initiate a pull request in Bitbucket. This is done by going 
into the Toolkit BitBucket project and clicking the pull request button -- 
make sure you select the correct source and destination branches. By default, 
the destination branch is set up to be the develop branch. So, in most cases, 
you won't have to do anything special. You must also select appropriate team
members to review changes. Our Bitbucket project is set up to require at least 
one other developer to approve the pull request.

.. important:: **You cannot approve your own pull request.**

When your pull request is approved (see :ref:`review-label` for more 
information), you merge your changes to the develop branch by clicking the 
"merge" button in Bitbucket. If there are no merge conflicts, the merge will 
proceed and you are done.

If there are conflicts, Bitbucket will not allow the merge to be done. 
You must resolve the conflicts first. The preferred way to do this is to go 
into your branch and do the following::

  $ git fetch origin
  $ git merge origin
  $ git commit
  $ git push

The 'fetch' command pulls changes from the remote branch into your local branch.
Running the 'merge' command will show which files have conflicts. Conflicts
are resolvers as described in the previous step. After all conflicts are
resolved, run the 'commit' and 'push' commands as usual. Lastly, complete
the merge in Bitbucket by clicking the merge button.

.. important:: **To keep things tidy, delete your topic branch in Bitbucket after it is merged into develop if you do not need it for further development. This can also be done by selecting the option in Bitbucket before doing the merge.**

Checking out an existing branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Any existing branch can be checked out from the Git repository and cloned from,
etc. Here are some useful commands::

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

