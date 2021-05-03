.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core concepts
******************************************************

Sidre provides five main classes: Datastore, Buffer, Group, View, and 
Attribute. In combination, these classes implement a data store with a tree 
structure to organize data in a hierarchy:

* DataStore is the main interface to access a data hierarchy.
* Buffer describes and holds data in memory.
* Group defines parent-child relationships in a hierarchical tree data structure and provides access to file I/O operations.
* View provides a virtual description of data and access to it.
* Attribute allows a program to attach metadata to View objects for processing data selectively. 

The following sections summarize the main interface features and functionality
of these Sidre classes.

.. note:: Interfaces for each of these classes are provided natively in C++, 
          C, and Fortran. 

.. toctree::
   :maxdepth: 2

   datastore
   buffer
   group
   view 
   attribute

