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
Software Development Tools and Usage
======================================================

We use a variety of tools for software development. Our tool philosophy has
three main tenets:

  * Employ robust, commonly-used tools. Don't re-invent something that already exists
  * Apply tools in ways that are easy for non-experts
  * Strive for automation and reproducibility

The main tools we use are listed below. Details about how we use 
them and helpful information about getting started are provided 
in the sections that follow.

* We use our `Confluence space <https://lc.llnl.gov/confluence/display/ASCT/ASC+Simulation+CS+Toolkit+Home>`_ for team discussion, planning, maintaining meeting notes, etc.
* Our Git repository contains Toolkit source code, build configurations, scripts, test suites, documentation, etc. The repository lives in our `Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.
* We use our `JIRA project <https://lc.llnl.gov/jira/projects/ATK>`_ for issue tracking.
* We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for continuous integration and automated testing.


--------------------------
Build, Link, Triumph (BLT)
--------------------------

Our build system, called *BLT*, is maintained in `its own repo <https://lc.llnl.gov/bitbucket/projects/ATK/repos/blt/browse>`_ in our 
Bitbucket project. BLT provides a "common sense" setup based on CMake for 
configuring and building the Toolkit code. It also enables software development 
tool integration via *make targets*. BLT has built-in support for the following
tools, all of which we use for Toolkit development:

  Documentation
    *Doxygen* (source code docs) and *Sphinx* (user docs)
  Unit testing
    *CTest* (test orchestration), *Google Test* (C/C++ unit tests), *Fruit* (Fortran unit tests)
  Code Health
    *Uncrustify* (code style), *gcov* and *lcov* (code coverage), and *Valgrind* (memory checking)
  Benchmarking
    *Google Benchmark*

See **BLT documentation (add link)** for more information.  

We use `Spack <https://github.com/LLNL/spack>`_ to manage and build the 
third-party libraries on which the Toolkit depends.

The Toolkit **Quick Start Guide (add link)** describes how to build the
code and third-party libraries.

.. note :: BLT is supported as a standalone product and used by other 
           software projects.


--------------------------------------
Git/Bitbucket
--------------------------------------

This section provides some information about getting started with Git and 
Bitbucket and describes operations related to topic branch development 
on the CS Toolkit project. Our Git repository lives in our 
`Bitbucket project <https://lc.llnl.gov/bitbucket/projects/ATK>`_.

If you are new to the Git or want to brush up on its features, there are 
several good source of information available on the web:

  * `Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of useful stuff.
  * The `Git Docs <https://git-scm.com/docs/>`_ is a complete reference for Git commands and options. It also provides soem *cheat sheets* you can download.
  * `Learn Git Branching <http://learngitbranching.js.org/>`_ is nice for visual, hands-on learners. 
  * The e-book `Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an excellent overview guide to using Git effectively.

To make Git easier to work with, you can define aliases in your shell
environment to do things like modify your prompt to show which git branch you
are on. If you are a csh/tcsh user, for example, you can add the following to
the file (e.g., .cshrc) that defines your profile::

   alias __git_current_branch 'git rev-parse --abbrev-ref HEAD >& /dev/null && echo "{`git rev-parse --abbrev-ref HEAD`}"'
   alias precmd 'set prompt="%n@%m>`__git_current_branch` "'

See also 
`Git scripts <https://github.com/git/git/tree/master/contrib/completion>`_ 
and elsewhere for useful scripts that folks have written to enable features
like tab-autocompletion for Git commands.

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

It is worth re-emphasizing a fundamental principle of the Gitflow 
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

  Step 1 -- Create a topic branch. 
    A topic branch name should include your 
    user id and a brief description indicating the purpose of the branch. We 
    typically label such branches using "feature", "bugfix", etc. to make it 
    clear what type of work is being performed on a branch. For example,::

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

  Step 2 -- Edit Files.
    After the topic branch is created, and you've pushed 
    it to Bitbucket, perform your development; i.e., edit files, add files, etc. 
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

  Step 3 -- Keep current with develop.
    If you will be working on your branch 
    for a while, it is a good idea to merge from the develop branch to your topic 
    branch regularly to avoid getting too far out of sync. Otherwise, you may have 
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
    
    Alternatively, you can use a tool to help resolve your conflicts. The 
    'git mergetool' command helps you run a merge tool. One such tool is the 
    "meld" tool, which is very powerful and intuitive. Diff tools like "tkdiff" 
    are also helpful for resolving merge conflicts.
    
    .. important:: **Git will not let you commit a file with merge conflicts.**


  Step 4 -- Create a pull request.
    When your work is complete, and you are 
    ready to merge your topic branch to the develop branch, you must initiate a 
    pull request in Bitbucket. This is done by going 
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
    proceed and you are done. If there are conflicts, Bitbucket will tell you
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
    'push' commands as usual::
    
      $ git commit
      $ git push
    
    Lastly, complete the merge in Bitbucket by clicking the merge button.
    
    .. important:: **To keep things tidy, please delete your topic branch in 
                   Bitbucket after it is merged if you no longer need it for 
                   further development. Bitbucket also provides a button to click  
                   on to do this after the merge is complete.**

Checking out an existing branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When working on multiple branches, or working on one with someone else on
the team, you will need to checkout a specific branch. Any existing branch 
can be checked out from the Git repository and cloned from, etc. Here are 
some useful commands::

  $ git fetch
  $ git branch -a
  $ git checkout <branch name>

The 'fetch' command retrieves new work committed by others on branches you may
have checked out, but *without merging* those changes into your local
copies of those branches. The 'branch' command lists all available remote 
branches. The 'checkout' command checks out
the specified branch into your local working space. 

.. note:: **You do not give the '-b' option when checking out an existing branch. 
          This option is only used when creating a new branch.**

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

  Open.
    Every issues starts out as an open issue. An open issue can 
    be assigned to someone or unassigned. When an issue is assigned, this 
    means that the assignee owns the issue and is responsible for working 
    on it. An open issue that is unassigned has not been been discussed or 
    reviewed, or no decision to act on it has been made. Typically, an open 
    issue means that it is not being worked on.
  In Progress.
    An issue in progress is one that is actively being worked on.
  Closed.
    When an issue is closed, work on it has been completed, or 
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
fields properly greatly helps other team members search through project issues
to find what they are looking for. Note that issue fields marked with a red 
asterisk are required. The others are not required, but may be used to include 
helpful information. The main issues we use regularly are:

  Project
    The CS Toolkit will show up as the default. You shouldn't need
    to change this.
  Issue Type
    We use only three issue types: *Bug*, *New Feature*, and
    *Task*. A bug is something broken that needs to be fixed. A new feature
    is something to add that increases functionality, enhances an interface,
    etc. Task is a "catch-all" issue type for any other issue.
  Summary
    Provide a short descriptive summary. A good (and brief)
    summary makes it easy to scan a list of issues to find one you are 
    looking for.
  Priority
    Select an appropriate issue priority to impart its level 
    of importance or urgency. Clicking on the question mark to the right of
    the priority field provides a description of each option.
  Components
    Each issue is labeled with the Toolkit component it 
    applies to. Other "component" labels indicate build system issues, 
    documentation issues, etc. 
  Assignee
    Unless you are certain which team member should be assigned
    the issue, choose 'Unassigned'. This will indicate that the issue requires
    discussion and review before it is assigned. The default assignee is the
    owner of the component you chose earlier if you make no choice.
  Reporter
    Unless you explicitly enter someone in this field, you, as
    the issue creator, will be the reporter. This is the correct choice in
    almost all cases.
  Description
    The description field should be used to include important
    details about the issue that will help the developer who will work on it.
  Environment
    The environment field can be useful when an issue affects a particular
    compiler or platform.

You may also use the other fields that appear if you think they will help
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


Issue assignee
^^^^^^^^^^^^^^^

Note that an assigned issue can be assigned to someone else to work on it.
An assigned issue can also be set back to 'Unassigned' if it needs further 
discussion by the team.

JIRA tips
^^^^^^^^^^

Here are some links to short videos (a couple of minutes each) that 
demonstrate how to use JIRA features:

   * `JIRA Instant Search Bar Demo <https://www.youtube.com/watch?v=ZmACxhzXLco&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=3>`_
   * `JIRA System Files Demo <https://www.youtube.com/watch?v=O08oySq043w&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=4>`_
   * `Creating and Editing JIRA Issues <https://www.youtube.com/watch?v=EsQ__dR6Nrw&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=5>`_


--------------------------------------
Bamboo Continuous Integration
--------------------------------------

We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for 
continuous integration and automated testing. We maintain a collection of
test plans for performing automated and manual builds, tests, and other
code health monitoring tasks.

Bamboo Agent Notes
^^^^^^^^^^^^^^^^^^^
The Bamboo server hands our scripts to it's associated 'agents' on the various clusters.
Each bamboo agent needs to be approved by an LC Atlassian admin in order to start executing Bamboo plans.  
The Atlassian admin will take care of associating your approved agent with your project and plan(s).

Restarting the Agent:
 On occasion, the agent can die.  This results in bamboo jobs being queued and stalled until the agent is restarted.  
 You must have access to the toolkit 'ATK' user to restart the agent. 
 To manually restart the CZ agent: ::

  $ ssh cab687 xsu atk
  $ cd /g/g16/atk/bambooAgent/asctoolkit.cab.llnl.gov
  $ ./bin/bamboo-agent.sh stop/start

.. note :: 
   Bamboo agents are created on and approved to run on specific nodes of a machine 
   and admin approval is required to create agents on a different node.
   The project has approved agents on the following nodes: 
       
     * CZ CHAOS: cab687
     * RZ CHAOS: rzalastor1
     * RZ TOSS 3: rzgenie2 
     * BGQ: vulcanlac3 

There are cron jobs on CZ and RZ that attempt to restart our agents every hour. 


You can view the cron jobs on the CZ using::

 $ ssh cab687 xsu atk
 $ cd /g/g16/atk/bamboo
 $ crontab -l czcrontab.txt 

And you can view the jobs on the RZ using::

 $ ssh rzalastor1 xsu atk
 $ cd /g/g16/atk/bamboo
 $ crontab -l rzcrontab.txt 


Quick setup for adding additional agents::


  $ pick a node, for example if we are to create a bamboo agent on rzgenie for asctoolkit
  $ atk@rzgenie2 ~/bamboo:/collab/usr/global/tools/bamboo/install-agent asctoolkit chang28@llnl.gov
  $ follow the instructions

.. note:: 
    After the agent is created, please contact Atlassian admin and get ready to start the bamboo agent. 
    Make sure you have a test plan set up to be attached to the agent. 



Agent Configuration:
  All of your Bamboo plan jobs are found in your build directory (all that are using the same agent, that is).  By default, this is under the directory where you started your agent.

To specify an alternative home directory, edit the wrapper.conf and restart your agent::

  $ vi <your-agent-home>/conf/wrapper.conf
  $ # change the following setting to the path you want your builds to run under
  $ wrapper.java.additional.1=-Dbamboo.home=/usr/workspace/wsrzc/atk/bamboo/asctoolkit-rzgenie2-1 (path to the build_dir
  $ restart the agent


Steps to Configure Bamboo Test Plan on a new system:

.. System could be a new architecture such as BGQ, or a new OS like TOSS3.
   I would describe the process that I used to set up BGQ test plan on bamboo.
   BGQ already has an agent in place on Vulcan.
..

  1. First we need a bamboo agent on the new system.  
  2. After the agent is up and running, we need to make sure the Third Party Libraries (TPL) are built. 
     TPL needs to happen before we can build the Asctoolkit code (for example, cmake needs to be ready). 
     To set up a new system, modify the ``compilers.yaml`` script under ``scripts/uberenv``. 
     A successful TPL build would generate host configuration files for each compiler defined in ``compilers.yaml``.
  3. The next step is to create a python script similar to ``llnl_cz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py``. 
     The script is located in ``scripts/uberenv/llnl_install_scripts``.

Automated plans
^^^^^^^^^^^^^^^^

Currently, we have the following test plans on CZ:

  Build and Test Develop Branch (all compilers, nightly)
    This is done on a nightly basis on the develop branch. 
  Build and Test Master Branch (weekly, all compilers)
    This is done on a weekly basis on the master branch. 
  Build and Test Branch (all compilers, on-demand)
    This is done manually on the master branch. 
  Develop and Test TPL (weekly)
    This is done on a weekly basis on the develop branch. 
  Smoke Test(On-Demand)
    This is done manually on the develop branch. 


Currentl,y we have the following test plans on RZ:

  Build and Test Develop Branch (all compilers, nightly, rzalastor)
    This is done on a nightly basis on the develop branch. 

Plan and Branches
^^^^^^^^^^^^^^^^^

To add a repository to a plan:

  1. Select Actions -> Configure Plan
  2. Select the Repositories tab
  3. Click the Add Repository button.
  
     Basic options:
       * Repository Host is "Bitbucket / Stash" (the cz server can also pull from Github)
       * Server is CZ Bitbucket (only option available)
       * Repository "ASC Simulation CS Toolkit / ASCToolkit"
       * Select the branch
     
     Advanced Options:
       * Default is to use shallow clones
       * Have to explicitly select 'Use Submodules', if you want them
       * Enable a quiet period to aggregate multiple commits before building
       * Can enable a quiet period or add a regular expression to exclude particular changesets
      
  4. Add a "Source Code Checkout" step to the plan's tasks to pull the latest code

To create plans that use the branches feature:

  ASC Toolkit has a nightly build plan that uses the develop branch as it's primary repository.  
  If you want to run the same plan on branches of this repository they can be set up a few different ways, 
  selecting specific branches in the repository and/or create branch plans for branches matching a regular expression.  
  The branches will then inherit all of the stages and jobs of the parent plan without having to duplicate the plan, 
  so any modifications to the parent plan apply to all associated branches.
  Link: `Atlassian article on feature branches <https://www.atlassian.com/continuous-delivery/building-feature-branches-with-bamboo>`_
    
  The quick way to enable branch plans:
  
  * Select Actions -> Configure Plan 
  * Select the Branches tab
  * Click the Create Plan Branch button in the Branches section (first section of the branches configuration)
  * Select the branches you would like to execute the parent plan on (this includes the triggers for the parent plan)
  * Select "Enable Branches" to make the branch plans active

Use a regular expression for your branch plan:

  * This would be useful to enable the sprint plans w/out having to worry about the sprint number
  * Also on the Branches tab of the plan configuration
  * In the Automatic branch management section
  * Select "Create plan branches for matching new branches" from the New Branches listbox.  
    * Add a regular expression in the 'Match name' text box (something like "/sprint\/([0-9]*)/" or "/feature\/")
    * Determine if you want Bamboo to delete plan branches after a period of time or a period of inactivity.  These are both set to do not delete by default, but once you selct  the "Create plan branches for matching new branches" option they are set to automatically delete.
    * Branch merging is disabled by default (this would automatically merge branches if tests are successful)
    * IRA feature branches is selected by default, so if you enable the branches on this page, Bamboo will automatically create plan branches for branches that contain a JIRA ticket in the name.
    * Select triggers - either inherit the parent plan triggers or run the branch plan manually.

To execute a test plan/branch from command line:
  * Use this python script to execute a test plan /branch from a command line: /usr/bin/python ./queue_build.py
  * Use Usetn key can be found in this directory: login vulcanlac3 as atk, go to /g/g16/atk/bambooWorkspace/asctoolkit.cab.llnl.gov/xml-data/build-dir
  * Plan key can also be found from the test plan execution log file.

Who Can do What
^^^^^^^^^^^^^^^^
Bamboo allows certain tasks to be down with an elevated privilege. If one does not have the privilege, he/she cannot even see the screen/button. That causes major confusion among users. This cheat sheet is intended to provide guide line of what tasks can only be done by Admin, and what tasks can be done by Admin and users alike.

Tasks that can only be done by Atlassian admin:

  * Delete a plan.
  * Delete a job of a plan
  * Configure branches
  * Approve New Bamboo agent
  * Assign agent to a plan


Tasks that can be done by everyone:

  * Create a plan.
  * Configure a plan
  * Limit the job to run on Agent
  * Review agent log,  located at /g/g16/atk/bambooAgent/asctoolkit.cab2.llnl.gov/atlassian-bamboo-agent.log (asctoolkit.cab.llnl.gov)

