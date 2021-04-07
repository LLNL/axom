.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _buffer-label:

==========
Buffer
==========

A Sidre Buffer object holds an array of data *described* by a data type and 
length. The data owned by a Buffer is unique to that Buffer object; i.e., 
Buffer objects do not share data. 

A Buffer can be created without a data description and then described 
later in a separate operation, or it can be described when it is created. 
In either case, data description and allocation are distinct operations. This
allows an application to create buffers it needs, then assess the types and
amount of data they will hold before deciding how and when to allocate data.

.. note:: * Buffer objects can only be created and destroyed using DataStore 
            methods. The Buffer constructor and destructor are private 
            (see :ref:`datastore-label`).
          * Each Buffer object has a unique integer identifier generated when it
            is created. If you want to interact with a Buffer object directly,
            you must keep a pointer to it or note its id so that you can 
            retrieve it from the DataStore when needed.

Buffer objects are used to hold data for Sidre View objects in most cases.
Each Buffer object maintains references to the Views that refer to its data. 
These references are created when a Buffer object is attached to a View, or
data is allocated through a View. Data stored in a Buffer may be accessed 
through a View object or through the Buffer directly. See :ref:`view-label` 
for more information about Views.

The Buffer interface includes the following operations:

 * Retrieve the unique id of the Buffer object.
 * Query whether a Buffer is *described* or *allocated*.
 * Describe Buffer data (type and number of elements).
 * Allocate, reallocate, deallocate Buffer data.
 * Copy a given number of bytes of data from a given pointer to a Buffer
   allocation.
 * Get data held by a Buffer as a pointer or conduit::Node::Value type.
 * Get information about data held by a Buffer: type, number of elements,
   total number of bytes, number of bytes per element, etc.
 * Retrieve the number of Views the Buffer is attached to.
 * Copy Buffer description and its data to/from a conduit::Node.


