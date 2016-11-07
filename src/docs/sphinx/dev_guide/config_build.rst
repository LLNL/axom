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

.. _configbuild-label:

======================================================
Configure and Build
======================================================

This section describes more advanced CS Toolkit build and configuration 
tasks that all developers should be aware of and be able to perform. 
More basic build instructions can be found in the User Quick Start Guide
**add link to that guide**. See :ref:`repoclone-label` for information 
about accessing the code.


.. _hostconfig-label:

------------------
Host-config Files
------------------

We use *host-config* files to track build configurations we support and 
maintain reproducibility. We maintain a collection of such files in the 
'host-configs' directory for platforms and compilers we support. 
When passed to CMake, using the '-C' option, a host-config file initializes 
the CMake cache with the configuration specified in the file. 

.. note :: Need to describe how the host config files get generated and how
           to generate new ones.


--------------------------
Make Targets
--------------------------

Anything new to add, or is quickstart guide sufficient?


.. _tpl-label:

--------------------------
Third-party Libraries
--------------------------

Describe how to run the scripts to install third-party libraries for 
testing different versions locally on a branch and for installing new
libraries for the team to use...

Building and installing TPLs for all compilers on LC CHAOS platforms (CZ)::

   $ python ./scripts/uberenv/llnl_install_scripts/llnl_cz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py

Questions we need to answer include:

  * How does one add a new compiler or platform to the mix?
  * How does one build a new set of TPLs with for a single platform or compiler
    for testing?
  * What is the procedure for changing versions of one or more TPLs?
  * How do we keep things straight when using different TPL versions for 
    different branches?
  * How to use the scripts for team TPL support vs. local development 
    experimentation?
  * Others?

.. note :: Pull in content from ../web/build_system/thirdparty_deps.rst ...
           fill in gaps and make sure it it up-to-date...
           

Our lc_install scripts have diverged from their original intended purpose of allowing anyone on the team to easily build TPLs on the LC platforms we care about. Recall, we want as close to a 1-button push build solution as we can feasibly support. To achieve this dream – we need scripts that don't take any arguments and just work when called. Of course things will break occasionally, the process itself should be as simple as possible.  Here is a strawman for an ideal process:

    Submit a batch job to build and check a set of tpl installs:
        build the TPLs for all compilers to a unique shared install directory 
        patch the generated host config files with extra settings from revision controlled "manual.edits.txt" files
        build and test the current toolkit checkout against all host-configs
    If you need to update the host-configs for a branch (develop, your work, etc) to use the new libs:
        Copy the generated host configs into the branch and follow the right procedure for that branch. 
            For example:
                (for develop, use branch and a PR to merge -- hopefully the PR will provide evidence from step 2 that the new TPLs work)
                (for your branch, other branches – commit, or use a PR, whatever is appropriate)
    Remove stale TPL sets (the group needs to coordinate on when this can happen)
        (we could write a tool to check if any branches are using host configs that link them to a given unique shared install directory)

To move towards this solution, this PR includes a simplified set of scripts that allow us to build TPLs for the all compilers we want to use on the chaos5 systems on the CZ and RZ.  This PR does not change the existing scripts.

On both the RZ and CZ, these scripts install to the new /usr/workspace file system. This file system has a much higher space quota than the /usr/gapps/ file system, however,  our team's shared toolkit directory exists in different places on the CZ and the RZ, so the PR includes separate scripts for the CZ and RZ.
Using a batch job to build and test chaos5 TPLs on the CZ and RZ

CZ: (from surface, cab, etc)
# from asctoolkit git root
cd scripts/uberenv/lc_install_scripts
msub msub_llnl_cz_chaos5_all_compilers.sh

This submits a job to the batch system that asks for 1 node for 8 hours to do the builds.

The batch job runs llnl_cz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py to build TPL sets for all compiler combos we want for chaos5 on the CZ.

This python script creates a time-stamped directory under /usr/workspace/wsa/toolkit/thirdparty_libs/builds/ and uses uberenv to install TPLs for a list of Spack specs.

It will also patch the generated host config files using "manual.edits.txt" files  and build the toolkit against all host configs

The batch job will log its output to a file created using the pattern: m.out.r.uberenv.chaos5.all.compilers.{job_id}.{hostname}.txt
RZ: (from rzmerl)
# from asctoolkit git root
cd scripts/uberenv/lc_install_scripts
msub msub_llnl_rz_chaos5_all_compilers.sh

This submits a job to the batch system that asks for 1 node for 4 hours to do the builds.

The batch job runs llnl_rz_uberenv_install_chaos_5_x86_64_ib_all_compilers.py to build TPL sets for all compiler combos we want for chaos5 on the RZ.

This python script creates a time-stamped directory under /usr/workspace/wsrzc/toolkit/thirdparty_libs/builds/ and uses uberenv to install TPLs for a list of Spack specs.

It will also patch the generated host config files using "manual.edits.txt" files  and build the toolkit against all host configs

The batch job will log its output to a file created using the pattern: m.out.r.uberenv.chaos5.all.compilers.{job_id}.{hostname}.txt
