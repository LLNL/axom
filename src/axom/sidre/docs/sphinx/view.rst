.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _view-label:

==========
View
==========

A Sidre View describes data and provides access to it. A View can describe
a (portion) of a data allocation in **any way that is compatible** with the
allocation. Specifically, the allocation must contain enough bytes to support
the description. In particular, the data type of a View description need
not match the types associated with the allocation employed by other Views
into that data.

.. note:: View objects can only be created and destroyed using Group methods 
          provided for this. The View constructor and destructor are private.

Each View object has a name and is owned by one Group in a Sidre Group 
hierarchy; its *owning* Group. A View maintains a pointer to the Group that 
owns it.

.. note:: * The name (string) of a View **must be unique** within its
            owning Group.
          * A View has a unique integer identifier within its owning group, 
            which is generated when the View is created.
          * Views in a Group can be accessed by name or integer id.

A View object can describe and provide access to data referenced by a 
pointer in one of four ways described below. In that case, a View data 
description includes: a data type, a length (number of elements), an offset 
and a stride (based on the pointer address and data type). 

  * A View can describe (a subset of) data owned by a pre-existing Buffer.
    In this case, the Buffer is manually *attached* to the View and the
    View's data description is applied to the Buffer data. Buffer data can be 
    (re)allocated or deallocated by the View if and only if it is the only 
    View attached to the Buffer. **In general, a Buffer can be attached
    to more than one View.**
  * A View description can be used to allocate data for View using a View
    ``allocate()`` method similar to Buffer data description and allocation 
    (see :ref:`buffer-label`). In this case, the View is usually exclusively 
    associated with a Buffer and no other View is allowed to (re)allocate or 
    deallocate the data held by the Buffer.
  * A View can **describe** data associated with a pointer to an *external* 
    data allocation. In this case, the View cannot (re)allocate or deallocate 
    the data. However, all other View operations can be applied to the data
    in essentially the same ways as the previous two cases.
  * A View can hold a pointer to an undescribed (*opaque*) data pointer. In 
    this case, the View knows nothing about the type or structure of the data; 
    it can only provide access to it. A user is entirely responsible for 
    casting the pointer to a proper type, knowing the size of the data, etc.

A View may also refer to a scalar quantity or a string. Such Views hold their
data differently than the pointer cases described above.

Before we describe the Sidre View interface, we present some View concepts
that describe various *states* a View can be in at any given time. Hopefully,
this will provide some useful context for the method descriptions that follow.

The key View concepts that users should be aware of are: 

  * View data description (data type, number of elements, stride, offset, etc.)
  * View data association (data lives in an attached Sidre Buffer object,
    accessed via external pointer, or is a scalar or string owned by the View)
  * Whether the View data description is applied to the data

The table below summarizes View data associations (rows) and View states with 
respect to that data (columns).

.. figure:: figs/sidre-view-states.png

   This table summarizes Sidre View *data associations* and *data states*. 
   Each row is a data association and each column refers to a data state.
   The True/False entries in the cells indicate return values of the
   View methods at the tops of the columns. The circumstances under which those
   values are returned are noted as well.

The three View data state methods at the tops of columns and their return 
values are:

  * **isDescribed()** returns true if a View has a data description, and
    false otherwise.
  * **isAllocated()** returns true if a View is associated with data, such as
    a non-null pointer, and false otherwise.
  * **isApplied()** returns true if the View has a data description and is
    associated with data that is compatible with that description, and the 
    description has been applied to the data; otherwise false is returned.

The rows indicate data associations; the View interface has methods to query
these as well; e.g., isEmpty(), hasBuffer(), etc. The associations are:

  * **EMPTY.** A View with no associated data; the View may or may not have
    a data description.
  * **BUFFER.** A View with an attached buffer; the View may or may not have 
    a data description and the Buffer may or may not be allocated and the
    description (if View has one) may or may not be applied to the Buffer data
    (if allocated).
  * **EXTERNAL.** A View has a non-null pointer to external data; the View
    may or may not have a data description and the description (if View has one)
    may or may not be applied to the external data.
  * **SCALAR.** View was created to hold a scalar value; such a View always
    has a valid data description, is allocated, and description is applied.
  * **STRING.** View was created to hold a string; such a View always
    has a valid data description, is allocated, and description is applied.

Note that there are specific consequences that follow from each particular
association/state that a View is in. For example, an EMPTY View cannot have an
attached Buffer. Neither can an EXTERNAL, SCALAR or STRING View. A View that
is EMPTY, BUFFER, SCALAR, or STRING cannot be EXTERNAL. Etc.

The following lists summarize the parts of the View interface:

.. note:: Most View methods return a pointer to the View object on which the
          method is called. This allows operations to be chained; e.g., ::

             View* view = ...;
             view->describe(...)->allocate(...)->apply(); 

.. _view-interface-label:

View Property Operations
-----------------------------

 * Retrieve the name or id of the View object.
 * Retrieve the View path name from the root of the tree or the path to the
   Group that owns it.
 * Get a pointer to the Group that owns the View.
 * Is View equivalent to another View; i.e., are names and data descriptions
   the same?
 * Rename a View.

Data Association Queries
--------------------------

 * Is View empty?
 * Does View have a Buffer attached?
 * Is View associated with external data?
 * Is it a scalar View?
 * Is it a string View?

Data State Queries
-------------------

 * Does View have a data description?
 * Is View data allocated?
 * Is View data description applied to data?
 * Is View opaque; i.e., it has an external pointer and no description?

Data Description Queries
--------------------------

 * Get type of data.
 * Get total number of bytes.
 * Get number of elements (total bytes / size of type).
 * Get number of bytes per data element (for type).
 * Get data offet.
 * Get data stride.
 * Get number of dimensions and shape of multi-dimensional data.
 * Get a conduit::Schema object that contains the View data description.

Data Management Operations
---------------------------

 * Allocate, reallocate, and deallocate View data.
 * Attach Buffer to View (with or without data description), 
   and detach Buffer from View.
 * Apply current View description to data or apply a new description.
 * Set View scalar value.
 * Set View string. 
 * Set external data pointer, with or without a data description. 

Data Access Methods
-----------------------

 * Get a pointer to View data, actual type or void*.
 * Get scalar value for a scalar View.
 * Retrieve pointer to Buffer attached to View.
 * Get a conduit::Node object that holds the View data.

Attribute Methods
-------------------

 * Query whether View has an Attribute with given id or name.
 * Get Attribute associated with a View by id or name.
 * Query whether Attribute has been set explicitly for View.
 * Reset Attribute with given id or name to its default value.
 * Set Attribute with given id or name to a given scalar value or string.
 * Retrieve scalar value or string of an Attribute.
 * Iterate over Attributes of a View.

I/O Operations
--------------
 
 * Copy View data description to a conduit::Node.
