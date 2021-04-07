.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _attribute-label:

==========
Attribute
==========

Sidre Attributes enable attaching metadata (strings and values) to Sidre 
Views to support queries (e.g., search for Views with a given attribute name),
outputting data for a subset of Views to files, and other ways an application
may need to selectively process Views in a Sidre DataStore hierarchy.

An Attribute is created with a string name and a default scalar or string value.
A default value can be changed later as needed.

.. note:: * Attribute objects can only be created and destroyed using DataStore
            methods. The Attribute constructor and destructor are private 
            (see :ref:`datastore-label`).
          * Each Attribute has a unique name and integer identifier. Either can
            be used to retrieve it from the DataStore.

Each Sidre View inherits all Attributes contained in the DataStore at their 
default strings or values. Then, an application may explicitly set any
Attribute on a View. The application may also query the value of a View 
Attribute, query whether the Attribute was explicitly set, or set the 
Attribute back to its default value. See :ref:`view-label`
for more information about Views.

The Attribute interface includes the following operations:

 * Retrieve the name and unique id of the Attribute object.
 * Set the scalar or string value of an Attribute.
 * Get the type of an Attribute's scalar value.

