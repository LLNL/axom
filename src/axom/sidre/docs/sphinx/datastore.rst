.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _datastore-label:

==========
DataStore
==========

A Sidre DataStore object provides the main access point for Sidre contents,
including the data managed by Sidre. In particular, a DataStore maintains the 
group at the root of the Sidre group hierarchy, a collection of
Buffer objects, and a collection of Attribute objects. Generally, the first 
thing a Sidre user does is create a DataStore; this operation also creates 
the root group. Apart from providing access to the root group, a DataStore 
object provides methods to interact with Buffer and Attribute objects. 

.. note:: Buffer and Attribute objects can only be created and destroyed 
          using DataStore methods noted below. Their constructors and
          destructors are private.

DataStore methods for Buffers support the following operations:

 * Create, destroy, and allocate data in Buffer objects
 * Query the number of Buffers that exist
 * Query whether a Buffer exists with given id
 * Retrieve Buffer with given id
 * Iterate over the set of Buffers in a DataStore

Please see :ref:`buffer-label` for more information about using Buffer objects.

DataStore methods for Attributes support the following operations:

 * Create and destroy Attributes
 * Query the number of Attributes that exist
 * Query whether an Attribute exists with given name or id
 * Retrieve Attribute with given name or id
 * Iterate over the set of Attributes in a DataStore

Please see :ref:`attribute-label` for more information about using Attribute
objects.
