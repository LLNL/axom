.. ##
.. ## Copyright (c) 2017-18, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

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

.. toctree::
   :maxdepth: 2

   datastore
   buffer
   group
   view 
   attribute

