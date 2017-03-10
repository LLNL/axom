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

****************************
Axom (a.k.a. The CS Toolkit)
****************************

.. note:: The project recently changed its name from "The CS Toolkit" to 
          "Axom". We have not completed converted everything on the project
          to the new name. The contents of this page, and linked pages, will
          continue to be updated as the conversion continues.

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

  *  Provide the CS infrastructure foundation of the ECP ATDM multi-physics application at LLNL
  *  Support current ASC and other production applications and as they continue to evolve
  *  Provide capabilities for LLNL research codes, proxy apps, etc. that simplify technology transfer from research efforts into production applications

**This web page is the main place to find information about Axom.**


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
  *  Lumberjack (Scalable parallel message logging and filtering)
  *  `Sidre: Simulation Data Repository <../../sidre_docs/html/index.html>`_
  *  `Spio (Sidre Parallel I/O) <../../spio_docs/html/index.html>`_
  *  Slam (Set-theoretic Lightweight API for Meshes)
  *  Quest (Querying on Surfaces Tool)
  *  Mint (Mesh data model)
  *  Primal (Computational geometry primitives)

--------------------------
Source Code Documentation
--------------------------

  *  `Axom <../../../doxygen/asc_toolkit_doxygen/html/index.html>`_
  *  `Common <../../../doxygen/common_doxygen/html/index.html>`_
  *  `Slic <../../../doxygen/slic_doxygen/html/index.html>`_
  *  `Lumberjack <../../../doxygen/lumberjack_doxygen/html/index.html>`_
  *  `Sidre <../../../doxygen/sidre_doxygen/html/index.html>`_
  *  `Spio <../../../doxygen/spio_doxygen/html/index.html>`_
  *  `Slam <../../../doxygen/slam_doxygen/html/index.html>`_
  *  `Quest <../../../doxygen/quest_doxygen/html/index.html>`_
  *  Mint
  *  Primal

Look for documentation to appear for new components as they are developed.


======================================================
Other Tools Application Developers May Find Useful
======================================================

Axom developers support other tools that can be used by software 
projects independent of the Axom. These include:

  *  `BLT <https://github.com/LLNL/blt>`_ (CMake-based buld system developed by the Axom team to simplify CMake usage and development tool integration)
  *  Shroud (Generation of native C and Fortran APIs from C++ code)
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

You can add or remove yourself from either of these lists via the 
`LLNL E-Mail List Manager (ListServ) <https://listserv.llnl.gov>`_


-------------- 
Chat Room
-------------- 

We also have a chat room on LLNL's Cisco Jabber instance called 
'Axom Dev'. It is open to anyone. You just have to log on to Jabber and
join the room.


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
Axom Release and Copyright License Information
======================================================

Fill this in when we have it..


.. toctree::
   :maxdepth: 3

