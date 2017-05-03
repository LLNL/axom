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

.. _bitbucket-label:

******************************************************
Git/Bitbucket: Version Control and Branch Development 
******************************************************

This section provides some information about getting started with Git and
Bitbucket and describes operations related to topic branch development
on the Axom project. Our Git repository lives in our
`Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.

If you are new to the Git or want to brush up on its features, there are
several good sources of information available on the web:

  * `Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of useful stuff.
  * The `Git Docs <https://git-scm.com/docs/>`_ is a complete reference for Git commands and options. It also provides some *cheat sheets* you can download.
  * `Learn Git Branching <http://learngitbranching.js.org/>`_ is nice for visual, hands-on learners.
  * The e-book `Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an excellent overview guide to using Git effectively.

=========
SSH Keys
=========

If you have not used Bitbucket before, you will need to
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_. This will make it easier for you to interact with our Git repository
without having to repeatedly enter login credentials. 

.. _repoclone-label:

=========================================
Cloning the Repo
=========================================

Before doing any work on the code, you must clone the repo into a local 
workspace. This is done by typing::

  $ git clone --recursive ssh://git@cz-bitbucket.llnl.gov:7999/atk/axom.git

.. note:: You don't need to remember the URL for the Axom repo above. It can be
          found by going to the Axom repo on our Bitbucket project and
          clicking on the 'Clone' action button that appears when you hover your
          mouse cursor over the ellipses on the top left of the web page.

The '--recursive' argument above is needed to pull the BLT build system, which
is a Git sub-module in Axom, into your local copy of the repo. In case you 
forget pass the '--recursive' argument to the 'git cone' command, you can
type the following commands after cloning::

  $ cd axom
  $ git submodule init
  $ git submodule update

Either way, the end result is the same and you are good to go.

After cloning, enter the top-level Axom directory and run the development
setup script we provide to ensure that your Git environment is configured
properly and client-side hooks we use are installed; i.e.,::

  $ cd axom
  $ ./scripts/setup-for-development.sh

You can also define aliases in your shell environment to do things like modify 
your prompt to show which git branch you are on. If you are a csh/tcsh user, 
for example, you can add the following to the file (e.g., .cshrc) that defines 
your profile::

   alias __git_current_branch 'git rev-parse --abbrev-ref HEAD >& /dev/null && echo "{`git rev-parse --abbrev-ref HEAD`}"'
   alias precmd 'set prompt="%n@%m>`__git_current_branch` "'

.. _topicdev-label:

=======================================
Topic Branch Development
=======================================

It is worth re-emphasizing a fundamental principle of the Gitflow
development model that we described in :ref:`gitflow-label`.

.. important:: **We never work directly on the develop or master branches.
               All development occurs on topic branches.**

When we refer to a *topic branch*, it could be a *feature branch*,
a *bugfix branch*, etc. The basic workflow for performing development
on a topic branch is:

  #. Create a topic branch off the develop branch and push the new branch
     to Bitbucket.
  #. Make changes and commit them to your branch in your local copy of the repo.
     Remember to push changes to the main repo on Bitbucket regularly for
     backup.
  #. If you are working on your topic branch for a while, it is a good idea
     to keep your topic branch current with the develop branch by merging 
     develop into your topic branch regularly.
  #. When your work is complete (including required testing, documentation, 
     etc.), create a pull request so others on the team
     can review your work. See :ref:`pullrequest-label`.

Here are some details about each of these steps.

Step 1 -- Create a topic branch
--------------------------------

    Most development occurs on a topic branch created off the develop branch.
    Occasions where a branch is created from another branch, such as a
    'hotfix' branch created off master, are described in :ref :`gitflow-label`.
    To create a branch in Git, you provide the '-b' argument to the 
    'git checkout' command. A topic branch name should include your user 
    id and a brief description indicating the purpose of the branch. 
    Typically, we label such branches using "feature", "bugfix", etc. to make 
    it clear what type of work is being performed on a branch. For example,::

      $ git checkout -b feature/<userid>/some_cool_new_feature
      $ git push -u

    You can also attach a JIRA issue number to the branch name. Then, Bitbucket
    will associate the issue with the commit when you merge your branch to the
    develop branch. For example,::

      $ git checkout -b bugfix/<userid>/jira-atk-<issue #>
      $ git push -u

    In each of these examples, the 'git push -u' command pushes the branch to
    the Bitbucket server and it will appear in the list of branches you and 
    other developers can see there.

Step 2 -- Edit files
--------------------------------

    After you've created a topic branch and pushed it to Bitbucket, perform 
    your development work on it; i.e., edit files, add files, etc. 
    Common commands you will use are::

      $ git add <file>
      $ git commit
      $ git push

    The 'add' command adds a file (or files) to be staged for commit. The 
    'commit' command commits staged files to your local copy of the repository.     The 'push' command pushes your commits to the topic branch in the main 
    Git repo. You could also do::

      $ git push origin

    This is equivalent to 'git push' if you specified the '-u' option when you
    originally pushed your topic branch you created it.

    Recall the Git environment setup script we recommended that you run after
    cloning the repo in the Section :ref:`repoclone-label` above. One of the
    Git pre-commit hooks that the script sets up applies formatting constraints
    on the commit message you provide when you execute the 'commit' command. The
    constraints are recommended Git practices that help make it easier to use
    various tools with the Git version control system. Specifically:

    * Commit message subject line is at most 50 characters
    * Subject line and main body of commit message are separated by a blank line
    * Main body of commit message is wrapped to 78 characters

Step 3 -- Keep current with develop
-------------------------------------

    If you will be working on your branch for a while, it is a good idea to 
    merge from the develop branch to your topic branch regularly to prevent 
    getting too far out of sync with the branch into which your work will be
    merged eventually. Otherwise, you may have many conflicts to resolve when 
    you are ready to merge your topic branch into develop and the merge could 
    be difficult.

    Before you begin the merge, make sure all outstanding changes to your topic
    branch are committed. Then, make sure your local repo is up-to-date with 
    the main develop branch by checking it out and pulling in the latest 
    changes; i.e.,::

      $ git checkout develop
      $ git pull

    Next, checkout your topic branch and merge changes in from the
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

    The section above the '=======' line are the file contents in the current 
    branch head (your topic branch). The lines below are the contents of the 
    develop branch that conflict with yours. To resolve the conflict, choose 
    the correct version of contents you want and delete the other lines.

    Alternatively, you can use a tool to help resolve your conflicts. The
    'git mergetool' command helps you run a merge tool. One such tool is the
    "meld" tool, which is very powerful and intuitive. Diff tools like "tkdiff"
    are also helpful for resolving merge conflicts.

    .. important:: **Git will not let you commit a file with merge conflicts.**

.. _createpr-label:

Step 4 -- Create a pull request
-------------------------------------

    When your work is complete, and you are ready to merge your topic branch 
    to the develop branch, you must initiate a pull request in Bitbucket. Go
    into the Axom Bitbucket project, select your branch, and click 
    `Create pull request` in the left column. Make sure you select the correct 
    destination branch. The default destination branch in our project is set 
    up to be the develop branch. So, in most cases, you won't have to do 
    anything special.

    You must also select appropriate team members to review changes. Our 
    Bitbucket project is set up to require at least one other developer to 
    approve the pull request.

    .. important:: **You cannot approve your own pull request.**

    When your pull request is approved (see :ref:`review-label` for more
    information), you merge your topic branch to the develop branch by 
    clicking the "merge" button in Bitbucket. If there are no merge conflicts, 
    the merge will proceed and you are done. If there are conflicts, Bitbucket 
    will indicate this and will not let you merge until all conflicts are
    resolved.

    The preferred way to resolve conflicts at this point is to go into your 
    topic branch and do the following::

      $ git fetch origin
      $ git merge origin

    The 'fetch' command pulls changes from the remote branch into your local
    branch. Running the 'merge' command will show which files have conflicts.
    Fix the conflicts as described in the previous step. After all conflicts 
    are resolved, run the 'commit' and 'push' commands as usual::

      $ git commit
      $ git push

    Lastly, complete the merge in Bitbucket by clicking the merge button.

    .. important:: **To keep things tidy, please delete your topic branch in
                   Bitbucket after it is merged if you no longer need it for
                   further development. Bitbucket provides an option to delete
                   the source branch of a merge after the merge is complete.**

================================
Checking Out an Existing Branch
================================

When working on multiple branches, or working on one with someone else on
the team, you will need to checkout a specific branch. Any existing branch
can be checked out from the Git repository and cloned from, etc. Here are
some useful commands::

  $ git fetch
  $ git branch -a
  $ git checkout <branch name>

The 'fetch' command retrieves new work committed by others on branches you may
have checked out, but *without merging* those changes into your local
copies of those branches. You will need to merge branches if you want changes
from one branch to be moved into another. The 'branch' command lists all 
available remote branches. The 'checkout' command checks out the specified 
branch into your local working space.

.. note:: **You do not give the '-b' option when checking out an existing 
          branch. This option is only used when creating a new branch.**

Here is a concrete example::

  $ git branch -a | grep homer
    remotes/origin/feature/homer/pick-up-bart
  $ git checkout feature/homer/pick-up-bart
    Branch feature/homer/pick-up-bart set up to track remote branch feature/homer/pick-up-bart
    Switched to a new branch 'feature/homer/pick-up-bart'

