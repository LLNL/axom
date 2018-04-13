******************************************************
An introductory example
******************************************************

As an introduction to the core concepts in Sidre and how they work, here is an
example where we construct the Sidre Datastore shown in the following figure:

.. figure:: figs/sidre_datastore_example.png
   :figwidth: 650px
   :alt: diagram of an example datastore

   ..

   +------------------------------------+-------------+
   | Symbol                             | Sidre class |
   +====================================+=============+
   | .. image:: figs/roundrectangle.png | Group       |
   |    :width: 120px                   |             |
   +------------------------------------+-------------+
   | .. image:: figs/rectangle.png      | View        |
   |    :width: 61px                    |             |
   +------------------------------------+-------------+
   | .. image:: figs/hexagon.png        | Attribute   |
   |    :width: 72px                    |             |
   +------------------------------------+-------------+

The diagram represents a Datastore, which contains  all Sidre objects and
provides the main interface to access those objects.
Rounded rectangles represent Sidre Group objects.  Each Group has a name,
zero or one parent Group, and zero or more child Groups (indicated by an
arrow from the parent to each child).  The Datastore provides exactly one
root Group (i.e. "/"): an application doesn't need to create the root.
The root Group is the only Group that has no parent.
Each Group also owns zero or more View objects, which are
shown as rectangles.  An arrow points from a Group to each View it owns.

A Sidre View object has a name and some data associated with it.
This example shows various types of data that
can be described by a View, including scalars, strings, and data in arrays (both
externally allocated and managed by Sidre Buffer objects).  Each array View
has a data pointer and describes data in terms of data type, number of elements,
offset, and stride.  Data pointers held by array Views are shown as
dashed arrows.

Datastores contain a collection of Buffer objects, shown as segmented
rectangles.

A Datastore also contains a list of Attributes.  Each Attribute is outlined
with a hexagon and defines a metadata label and
the default value associated with that label.  In this example, the Datastore
has Attributes "vis" (with default value 0) and "restart" (with default
value 1).  Default Attributes apply to all Views unless explicitly set for
individual Views.  In this example, the Views "temp" and "rho" have the
Attribute "vis" set to 1.

Other aspects of Sidre class use
are clarified in the C++ code shown next.  Sidre provides full C and Fortran
APIs that can also be used to generate the same result.

First, we create a tree of Group objects.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_creategroups_start
   :end-before: _first_example_creategroups_end
   :language: C++

The ``Group::createViewScalar()`` method lets an application store scalar values
associated with a Group.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_state_start
   :end-before: _first_example_state_end
   :language: C++

This example stores (x, y, z) node position data in one array of triples.
The array is managed through a Buffer object and three Views point into it.
C++ Sidre operations that create Buffers, Groups, and Views, as shown in the
following code, return a pointer to the
object that is created. This allows chaining operations.  (Chaining is supported
in the C++ API but not in C or Fortran.)

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_nodes_start
   :end-before: _first_example_nodes_end
   :language: C++

The next snippet shows how to set an Attribute on a View and how to create a View
that holds an external pointer described with data type, element count, offset, and stride.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_fields_start
   :end-before: _first_example_fields_end
   :language: C++

Here is how to retrieve data items out of the hierarchy.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _first_example_access_start
   :end-before: _first_example_access_end
   :language: C++

In the last section, the code accesses the arrays associated with the views "y",
"temp", and "region".  While "temp" and "region" have the default offset (0) and
stride (1), "y" has offset 1 and stride 3.  The pointer returned by
``View::getPointer()`` always points to the first element described by the
View (that is, the View takes care of the offset), but use of a stride other than 1
must be done by the code itself.

Unix-like path syntax using the slash ("/") delimiter is supported for
traversing Sidre Group and View hierarchies.  This usage is shown in the last
line of the code example above.  The getView() method call retrieves the View
named "region" in the Group "ext" that is a child of the "fields" Group.
Character sequences before the first slash and between two consecutive slashes
are Group names (describing parent-child relationships).  For this method, and
others dealing with Views, the sequence following the last slash is the name of
a View.  Similar path syntax can be used to retrieve Groups, create Groups and
Views, and so forth.
