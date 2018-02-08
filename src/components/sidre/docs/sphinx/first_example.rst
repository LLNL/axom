******************************************************
An introductory example
******************************************************

As an introduction to the core concepts in Sidre and how they work, here is an
example where we construct the Sidre Datastore shown in the following figure:

.. image:: sidre_datastore_example.png

Here, ovals represent Sidre Group objects.  Each Group has a name, 
zero or one parent Group, and zero or more child Groups (indicated by
arrows between ovals).  Note that the root Group (i.e. "/") is the only Group
that has no parent.  Each Group also owns zero or more View objects.  Views are
shown as rectangles and an arrow points from the owning Group to each View.
Datastores contain a collection of Buffer objects, shown as segmented
rectangles.  View references to Buffers or to external pointers are shown as
dashed arrows.

A Sidre View object has a name and some data
associated with it.  Here, we show various types of data that can be described
by a View, including scalars, strings, and arrays (both "external" and managed
by Sidre Buffer objects).  In the array Views ("x", "y", "z", "temp", "rho",
"region"), we include a triple indicating the number of elements, offset, and
stride for the View.  Other aspects of View descriptions are clarified in the
C++ code shown next.  Sidre provides full C and Fortran APIs that can also be
used to generate the same result.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: first_example_create_start
   :end-before: first_example_create_end
   :language: C++

Sidre operations that create Buffers, Groups, and Views return a pointer to the
object that is created. This allows chaining operations.  From the previous 
example:

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: first_example_chain_1
   :end-before: first_example_chain_2
   :language: C++

Many other Sidre operations return a pointer to the object on which the method
is called, which enables similar operation chaining.

Lastly, we show a few instances of how to retrieve data items out of 
the hierarchy.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: first_example_access_start
   :end-before: first_example_access_end
   :language: C++
  
Unix-like path syntax using the slash ("/") delimiter is supported for
traversing Sidre Group and View hierarchies.  This usage is shown in the last
line of the code example above.  The getView() method call retrieves the View
named "region" in the Group "ext" that is a child of the "fields" Group.
Character sequences before the first slash and between two consecutive slashes
are Group names (describing parent-child relationships).  For this method, and
others dealing with Views, the sequence following the last slash is the name of
a View.  Similar path syntax can be used to retrieve Groups, create Groups and
Views, and so forth.
