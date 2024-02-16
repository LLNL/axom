.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _misctasks-label:

********************************
Miscellaneous Development Items
********************************

This section describes various development tasks that need to be 
performed at times and which are not covered in other sections.


===================
Web Documentation
===================

Axom web-based documentation is hosted on our 
`Read the Docs project <https://readthedocs.org/projects/axom/>`_. 
Multiple versions are visible there, including the latest content on the 
develop branch (*latest*) and the main branch (*main*). The documentation 
that appears is automatically re-generated each time a change is pushed to 
a branch in the GitHub repository that is enabled to be displayed on the 
Read the Docs project. If you are modifying Axom documentation, you can enable 
the branch you are working on so that you can see what it looks like as you 
push changes to the branch. If your documentation changes are part of a GitHub
pull request, it is a good idea to enable the documentation for that branch
and put a link to it in the pull request summary. This makes it easy for 
reviewers to check over your changes.

.. note :: When you no longer need the documentation of your branch to be
           visible on Read the Docs (e.g., your pull request is merged), 
           please disable that branch on Read the Docs.


===================
Code Health Tools
===================

This section describes how to run code health tools we use.


Code Coverage
---------------

Setting up and running code coverage analysis...


Static Analysis
---------------

Setting up and running static analysis tools....


Memory Checking
----------------

Setting up and running memory checking tools....
