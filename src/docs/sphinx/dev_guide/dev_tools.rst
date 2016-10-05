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

  * Adopt commonly-used tools (don't invent our own if we don't have to)
  * Make the tools easy to use for non-expers
  * Focus on automation and reporoducibility

Our key tools are listed here. Details about how we use them are provided in 
the following sections.

The main interaction hub for developers is the Atlassian tool suite on the 
Livermore Computing Collaboration Zone (CZ). Access to our Atlassian project
spaces requires membership in the LC groups 'toolkit' and 'toolkitd'. The 
following list contains links to each of these spaces: 

* Our Git repository houses the Toolkit source code, build configurations, test suites, documentation, etc. The repository lives in a `Bitbucket project <https://https://lc.llnl.gov/bitbucket/projects/ATK>`_
* We use our `Confluence project space <https://lc.llnl.gov/confluence/display/ASCT/ASC+Simulation+CS+Toolkit+Home>`_ for team discussion, planning, maintaining meeting notes, etc.
* We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for issue tracking.
* We use our `Bamboo project <https://lc.llnl.gov/bamboo/browse/ASC>`_ for continuous integration and automated testing.

Our build system, called *BLT*, is maintained in its own repo in our 
Bitbucket project. **Add link to BLT documentation** BLT offers 
a "common sense" setup based on CMake. It manages our build environment 
(compilers, programming models - OpenMP, MPI, CUDA, etc., and third-party 
library locations) as well as our software development tool integration
via *makefile targets*. BLT has built-in support for:

* Documentation - *Doxygen* (source code docs) and *Sphinx* (user docs)
* Unit testing - *CTest* (test orchestration), *Google Test* (C/C++ unit tests), *Fruit* (Fortran unit tests)
* Code Health - *Uncrustify* (code style), *gcov* and *lcov* (code coverage), and *Valgrind* (memory checking)
* Benchmarking - *Google Benchmark*

BLT is supported as a standalone product and used by other software projects. 

We use `Spack <https:://github.com/LLNL/spack>`_ to manage and build the 
third-party libraries on which the Toolkit depends and a custom python 
script to bootstrap Spack.


--------------------------------------
Git/Bitbucket Basics
--------------------------------------

This section discusses basic Git and Bitbucket operations related to feature 
branch development on the CS Toolkit project. 

If you are new to the Git or want to brush up on some of its features, the 
`Atlassian Git Tutorial <https://www.atlassian.com/git/>`_ has a lot of good 
information, and `Learn Git Branching <http://learngitbranching.js.org/>`_ 
is nice for visual, hands-on learners. Also, the book 
`Pro Git, by Scott Chacon <https://git-scm.com/book/en/v2>`_ is an
excellent comprehensive guide to Git. To make Git easier to work with, folks 
have written some useful scripts. For example, see `Git scripts <https://github.com/git/git/tree/master/contrib/completion>`_ for scripts that enable 
tab-autocompletion for Git commands and set up you prompt to show which branch
you are on, etc.

If you've not used Bitbucket before, you will need to 
`create an SSH key <https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html>`_ and `add the key to your Bitbucket profile <https://confluence.atlassian.com/bitbucketserver/ssh-user-keys-for-personal-use-776639793.html>`_.
