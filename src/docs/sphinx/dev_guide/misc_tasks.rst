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

.. _misctasks-label:

======================================================
Miscellaneous Development Tasks
======================================================

This section describes various development tasks that need to be 
performed that are not covered in earlier sections.


------------------
Web Documentation
------------------

Build and install web documentation...

Shared LC web content location asctoolkit/src/docs/sphinx/web


------------------------------------
Third-party Library Installation
------------------------------------

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
