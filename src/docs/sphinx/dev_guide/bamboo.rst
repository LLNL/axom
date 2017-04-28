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

.. _citesting-label:

======================================================
Continuous Integration (Bamboo)
======================================================

We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for
continuous integration and automated testing. We maintain a collection of
test plans for performing automated and manual builds, tests, and other
code health monitoring tasks.


.. note:: This section needs work and cleanup....

Bamboo Agent Notes
^^^^^^^^^^^^^^^^^^^
The Bamboo server hands our scripts to it's associated 'agents' on the various clusters.
Each bamboo agent needs to be approved by an LC Atlassian admin in order to start executing Bamboo plans.
The Atlassian admin will take care of associating your approved agent with your project and plan(s).

Restarting the Agent:
 On occasion, the agent can die.  This results in bamboo jobs being queued and stalled until the agent is restarted.
 You must have access to the Axom 'atk' shared user account to restart the agent.
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


  $ pick a node, for example if we are to create a bamboo agent on rzgenie for Axom
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
     TPL needs to happen before we can build the Axom code (for example, cmake needs to be ready).
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


Currently, we have the following test plans on RZ:

  Build and Test Develop Branch (all compilers, nightly, rzalastor)
    This is done on a nightly basis on the develop branch.

Plan and Branches
^^^^^^^^^^^^^^^^^

To add a repository to a plan:

  1. Select Actions -> Configure Plan
  2. Select the Repositories tab
  3. Click the Add Repository button.

     Basic options:
       * Repository Host is "Bitbucket Server / Stash" (the cz server can also pull from Github)
       * Server is CZ Bitbucket (only option available)
       * Repository "Axom"
       * Select the branch

     Advanced Options:
       * Default is to use shallow clones
       * Have to explicitly select 'Use Submodules', if you want them
       * Enable a quiet period to aggregate multiple commits before building
       * Can enable a quiet period or add a regular expression to exclude particular changesets

  4. Add a "Source Code Checkout" step to the plan's tasks to pull the latest code

To create plans that use the branches feature:

  Axom has a nightly build plan that uses the develop branch as its primary repository.
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
  * Select "Create plan branches for matching new branches" from the New Branches list
    * Add a regular expression in the 'Match name' text box (something like "/sprint\/([0-9]*)/" or "/feature\/")
    * Determine if you want Bamboo to delete plan branches after a period of time or a period of inactivity.  These are both set to do not delete by default, but once you select  the "Create plan branches for matching new branches" option they are set to automatically delete.
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


