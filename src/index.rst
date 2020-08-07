.. ## Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

****
Axom
****

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

==========
Components
==========

   * Inlet: Input deck parsing and information storage/retrieval
   * Lumberjack: Scalable parallel message logging and filtering
   * Mint: Mesh data model
   * Primal: Computational geometry primitives
   * Quest: Querying on surface tool
   * Sidre: Simulation data repository
   * Slam: Set-theoretic lightweight API for meshes
   * Slic: Simple Logging Interface Code
   * Spin: Spatial indexes

=============
Documentation
=============

  * :doc:`Quickstart Guide <docs/sphinx/quickstart_guide/index>`
  *  `Source documentation <doxygen/html/index.html>`_

.. list-table::
   :align: center

   * - Core
     -
     - `Source documentation <doxygen/html/coretop.html>`_
   * - Inlet
     - :doc:`User Guide <axom/inlet/docs/sphinx/index>`
     - `Source documentation <doxygen/html/inlettop.html>`_
   * - Lumberjack
     - :doc:`User Guide <axom/lumberjack/docs/sphinx/index>`
     - `Source documentation <doxygen/html/lumberjacktop.html>`_
   * - Mint
     - :doc:`User Guide <axom/mint/docs/sphinx/index>`
     - `Source documentation <doxygen/html/minttop.html>`_
   * - Primal
     - :doc:`User Guide <axom/primal/docs/sphinx/index>`
     - `Source documentation <doxygen/html/primaltop.html>`_
   * - Quest
     - :doc:`User Guide <axom/quest/docs/sphinx/index>`
     - `Source documentation <doxygen/html/questtop.html>`_
   * - Sidre
     - :doc:`User Guide <axom/sidre/docs/sphinx/index>`
     - `Source documentation <doxygen/html/sidretop.html>`_
   * - Slam
     - :doc:`User Guide <axom/slam/docs/sphinx/index>`
     - `Source documentation <doxygen/html/slamtop.html>`_
   * - Slic
     - :doc:`User Guide <axom/slic/docs/sphinx/index>`
     - `Source documentation <doxygen/html/slictop.html>`_
   * - Spin
     - :doc:`User Guide <axom/spin/docs/sphinx/index>`
     - `Source documentation <doxygen/html/spintop.html>`_


============================
Component Level Dependencies
============================

Dependencies between components are as follows:

- Core has no dependencies, and the other components depend on Core
- Slic optionally depends on Lumberjack
- Slam, Spin, Primal, Mint, Quest, and Sidre depend on Slic
- Mint optionally depends on Sidre
- Quest depends on Slam, Spin, Primal, and Mint
- Inlet depends on Sidre, and Slic

The figure below summarizes the dependencies between the components.  Solid links
indicate hard dependencies; dashed links indicate optional dependencies.

.. graphviz:: docs/dependencies.dot


======================================================
Other Tools Application Developers May Find Useful
======================================================

Axom developers support other tools that can be used by software 
projects independent of the Axom. These include:

  *  `BLT <https://github.com/LLNL/blt>`_: CMake-based build system developed by the Axom team to simplify CMake usage and development tool integration
  *  `Shroud <https://github.com/LLNL/shroud>`_: Generator for native C and Fortran APIs from C++ code
  *  `Conduit <https://lc.llnl.gov/confluence/display/CON/Conduit+Home>`_: Library for describing and managing in-memory data structures

===================
Developer Resources
===================

  * :doc:`Developer Guide <docs/sphinx/dev_guide/index>`
  * :doc:`Coding Guide <docs/sphinx/coding_guide/index>`

================================
Communicating with the Axom Team
================================

--------------
Mailing Lists
--------------

The project maintains two email lists: 

  * 'axom-users@llnl.gov' is how Axom users can contact developers for questions, report issues, etc. 
  * 'axom-dev@llnl.gov' is for communication among team members. 


---------
Chat Room
---------

We also have a chat room on LLNL's Microsoft Teams called 'Axom'. They open to
anyone on the LLNL network. Just log onto Teams and join the room.


---------------------------------
Git repository and Issue Tracking
---------------------------------

The main interaction hub for Axom is on `Github <https://github.com/LLNL/axom>`_



======================================================
Axom Copyright and License Information
======================================================

Please see the :ref:`axom-license`.

Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.

LLNL-CODE-741217


.. toctree::
   :hidden:
   :maxdepth: 1

   docs/sphinx/quickstart_guide/index

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Component User Guides

   Slic (Simple Logging Interface Code) <axom/slic/docs/sphinx/index>
   Lumberjack (Scalable parallel message logging and filtering) <axom/lumberjack/docs/sphinx/index>
   Sidre (Simulation data repository) <axom/sidre/docs/sphinx/index>
   Slam (Set-theoretic lightweight API for meshes) <axom/slam/docs/sphinx/index>
   Spin (Spatial indexes) <axom/spin/docs/sphinx/index>
   Quest (Querying on surface tool) <axom/quest/docs/sphinx/index>
   Mint (Mesh data model) <axom/mint/docs/sphinx/index>
   Primal (Computational geometry primitives) <axom/primal/docs/sphinx/index>
   Inlet (Input decks) <axom/inlet/docs/sphinx/index>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Developer Resources

   docs/sphinx/dev_guide/index
   docs/sphinx/coding_guide/index
   docs/licenses
