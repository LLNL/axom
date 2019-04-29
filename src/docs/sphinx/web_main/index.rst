.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

****
Axom
****

.. note:: This page is under heavy construction.

Axom is a project in WCI/WSC that is funded by ECP/ATDM.
Its principal goal is to provide a collection of robust and flexible software 
components that serve as building blocks for LLNL simulation tools. The 
emphasis is on sharing core infrastructure software amongst applications 
rather than having different codes develop and maintain similar capabilities.

A key objective of Axom is to facilitate integration of novel, 
forward-looking computer science capabilities into LLNL simulation codes. 
Thus, a central function of Axom is to enable and simplify data exchange 
between applications and tools that Axom provides. To meet these 
objectives, developers of Axom components emphasize the following features 
in software design and implementation:

  * Flexibility to meet the needs of a diverse set of applications
  * High-quality, with well designed APIs, good documentation, tested well, high performance, etc.
  * Consistency in software engineering practices
  * Integrability so that components work well together and are easily adopted by applications

The main drivers of the Axom project are to:

  *  Provide the CS infrastructure foundation of the ECP/ATDM multi-physics application at LLNL
  *  Support current ASC and other production applications and as they continue to evolve
  *  Provide capabilities for LLNL research codes, proxy apps, etc. that simplify technology
     transfer from research efforts into production applications


=====================
Quickstart Guide
=====================

The `Axom Quickstart Guide <../../quickstart_guide_docs/html/index.html>`_ 
contains information about accessing the code, configuring and building, 
linking with an application, etc.


================================
Axom Software Documentation
================================

The following lists contain links to user guides and source code documentation
for Axom software components:

----------------------
Component User Guides
----------------------

  *  Slic (Simple Logging Interface Code for integrated applications)
  *  `Lumberjack (Scalable parallel message logging and filtering) <../../lumberjack_docs/html/index.html>`_
  *  `Sidre (Simulation data repository) <../../sidre_docs/html/index.html>`_
  *  `Slam (Set-theoretic lightweight API for meshes) <../../slam_docs/html/index.html>`_
  *  `Quest (Querying on surface tool) <../../quest_docs/html/index.html>`_
  *  `Mint (Mesh data model) <../../mint_docs/html/index.html>`_
  *  `Primal (Computational geometry primitives) <../../primal_docs/html/index.html>`_

--------------------------
Source Code Documentation
--------------------------

  *  `Axom <../../../doxygen/axom_doxygen/html/index.html>`_
  *  `Core <../../../doxygen/axom_doxygen/html/coretop.html>`_
  *  `Lumberjack <../../../doxygen/axom_doxygen/html/lumberjacktop.html>`_
  *  `Mint <../../../doxygen/axom_doxygen/html/minttop.html>`_
  *  `Primal <../../../doxygen/axom_doxygen/html/primaltop.html>`_
  *  `Quest <../../../doxygen/axom_doxygen/html/questtop.html>`_
  *  `Sidre <../../../doxygen/axom_doxygen/html/sidretop.html>`_
  *  `Slic <../../../doxygen/axom_doxygen/html/slictop.html>`_
  *  `Slam <../../../doxygen/axom_doxygen/html/slamtop.html>`_

Look for documentation to appear for new components as they are developed.

Dependencies between modules are as follows:

- Core has no dependencies, and the other modules depend on Core
- Slic optionally depends on Lumberjack
- Slam, Primal, Mint, Quest, and Sidre depend on Slic
- Mint optionally depends on Sidre
- Quest depends on Slam, Primal, and Mint

The figure below summarizes the dependencies between the modules.  Solid links
indicate hard dependencies; dashed links indicate optional dependencies.

.. graphviz:: dependencies.dot


======================================================
Other Tools Application Developers May Find Useful
======================================================

Axom developers support other tools that can be used by software 
projects independent of the Axom. These include:

  *  `BLT <https://github.com/LLNL/blt>`_ (CMake-based build system developed by the Axom team to simplify CMake usage and development tool integration)
  *  `Shroud <https://github.com/LLNL/shroud>`_ (Generator for native C and Fortran APIs from C++ code)
  *  `Conduit <https://lc.llnl.gov/confluence/display/CON/Conduit+Home>`_ (Library for describing and managing in-memory data structures) 


================================================
Resources for Axom Developers/Contributors:
================================================

  * `Axom Developer Guide  <../../dev_guide_docs/html/index.html>`_
  * `Axom Coding Guidelines  <../../coding_guide_docs/html/index.html>`_
  * `Axom Testing Coverage <https://lc.llnl.gov/toolkit/coverage/index.html>`_



======================================= 
Communicating with the Axom Team
=======================================

--------------
Mailing Lists
--------------

The project maintains two email lists: 

  * 'axom-users@llnl.gov' is how Axom users can contact developers for questions, report issues, etc. 
  * 'axom-dev@llnl.gov' is for communication among team members. 


-------------- 
Chat Room
-------------- 

We also have a chat room on LLNL's Cisco Jabber instance called 
'Axom Dev'. It is open to anyone on the llnl network. 
Just log onto Jabber and join the room.


-----------------
Atlassian Tools
-----------------

The main interaction hub for the Axom software is the Atlassian tool suite 
on the Livermore Computing Collaboration Zone (CZ). These tools can be 
accessed through the `MyLC Portal <https://lc.llnl.gov>`_.

Direct links to the Axom Atlassian projects/spaces are:

  * `Bitbucket project/git repository <https://lc.llnl.gov/bitbucket/projects/ATK>`_
  * `Jira issue tracker <https://lc.llnl.gov/jira/projects/ATK>`_
  * `Bamboo continuous integration <https://lc.llnl.gov/bamboo/browse/ASC>`_
  * `Confluence (primarily for developers) <https://lc.llnl.gov/confluence/display/ASCT>`_


--------------------
LC Groups
--------------------

Access to Axom projects/spaces in these Atlassian tools requires
membership in the `axom` group on LC systems. Please contact the team for
group access by sending an email request to 'axom-dev@llnl.gov'.


======================================================
Axom Copyright and License Information
======================================================

Please see the :ref:`axom-license`.

Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

LLNL-CODE-741217


.. toctree::
   :maxdepth: 3

   licenses
