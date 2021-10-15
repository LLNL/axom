.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _dataconcepts-label:

===========================
Data vs. Metadata Concepts
===========================

This section describes examples that illustrate some common Sidre usage 
patterns involving buffers, groups, and views. In particular, it shows
how data associated with a view is related to a buffer, how the lifetimes
of data associated with a view and the buffer owning the data are related,
and how the description of data in a view (metadata) is independent of the
the allocation state of the data in a view and a corresponding buffer. 

.. note:: Object and data creation and destruction is very explicit in Sidre
          to allow maximum flexibility to compose complex operations from 
          simpler ones. This should become apparent by reading the following 
          examples. 

The code snippets shown and described here exist in the file ``axom/src/axom/sidre/examples/sidre_data_vs_metadata.cpp``, which can be built and run to experiment with if you wish.

The starting point for the first three examples is a simple Sidre datastore 
group hierarchy in which the root group has two child groups named "A" and "B". 
This structure is generated as follows:

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex_datastore_initial_start
   :end-before: _ex_datastore_initial_end
   :language: C++

At this point, the datastore has zero buffers and the groups have no views
as we can see by the output of the code::

  Datastore start state....
        Num buffers in datastore: 0
        Num views in group A: 0
        Num views in group B: 0

--------------------------------
Example 1: One View, One Buffer 
--------------------------------

The first example illustrates a common Sidre usage pattern in which
a view is created in a group and the data associated with the view is managed
via the view. This case represents a one-to-one relationship between a view and
buffer.

We begin by creating a view "aview" in group "A" describing an integer array 
of length 10 and allocating the array, all in one method call. This creates a 
buffer in the datastore which holds the array data. Then, we get a pointer to 
the start of the array from the view and initialize the array values. To give 
some insight into how Sidre orks, we access and print various pieces of 
information about the state of the group, view, and buffer. We also print the 
value of the buffer array at slot 5 to confirm that the buffer is indeed 
holding the view's data. 

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_create_start
   :end-before: _ex1_oneview_onebuffer_create_end
   :language: C++

The output printed by the code is::

  After A_grp->createViewAndAllocate() call
        Num views in group A: 1
        Num elements in view: 10
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10

  After initialization of view array
        Value of elt 5 in buffer array (expect 7): 7

Next, we deallocate the view and show that its description remains intact; 
for example, the number of elements is still 10. The view is no longer 
allocated, but it is still attached to its buffer. We confirm this and that 
the buffer is no longer allocated as expected. 

Then, we allocate the view again. Since the buffer is still attached to the
view, a new buffer is not created in the datastore and the existing buffer is 
re-allocated. We verify this by showing that datastore still has one buffer.
The buffer and view are both shown to be allocated and the number of elements
in the view remains 10 since the view description has not changed.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_deallocalloc_start
   :end-before: _ex1_oneview_onebuffer_deallocalloc_end
   :language: C++

The output of the code is::

  After view deallocate call, the data no longer exists,
  but the view description remains.
        Num elements in view: 10
        View has buffer? 1
        Is view allocated? 0
        Num buffers in datastore: 1
        Is buffer allocated? 0
        Num views attached to buffer: 1
        Num elements in buffer array: 10

  After allocating view again...
        Num buffers in datastore: 1
        Is buffer allocated? 1
        Is view allocated? 1
        Num elements in view: 10 

Lastly, we destroy the view and its data with a single method call and verify
that the view and associated buffer no longer exist.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_destroy_start
   :end-before: _ex1_oneview_onebuffer_destroy_end
   :language: C++

The out of the code is::

  After destroyViewAndData() call
        Num views in group A: 0
        Num buffers in datastore: 0

---------------------------------
Example 2: Two Views, One Buffer 
---------------------------------

The second example illustrates a Sidre usage pattern in which multiple views 
are created to describe portions of data held in a shared buffer. This case 
represents a many-to-one relationship between views and a buffer. Before we 
start, we verify that the datastore contains no buffers and that the groups 
we created earlier contain no views.

We start by creating and allocating a buffer holding an array of doubles of 
length 10. We initialize the array so that each element has a value matching 
its position in the array; i.e., the values 0 through 9. Then, we create two 
views in group "A" each *attached* to the buffer. Next, we apply a data 
description containing an offset and stride to each view so that one view 
is associated with the even values in the buffer and the other is associated
with the odd values. Accessing the data pointer in each view and printing
the values shows that this is indeed the case.

We call a method to destroy the first view and its data, similar to the
last part of the first example. The view is destroyed. However, since the 
buffer that held its data is shared by the other view, the buffer and its
data remain intact. In addition, the data associated with the remaining view
is also untouched.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex2_twoviews_onebuffer_start
   :end-before: _ex2_twoviews_onebuffer_end
   :language: C++

The output of the code is::

  Datastore start state
        Num buffers in datastore: 0
       	Num views in group A: 0
        Num views in group B: 0

  After buffer allocation and attaching to views
       	Num buffers in datastore: 1
       	Buffer num elements: 10

       	aview1 data has even values:
       	aview1 num elements: 5
       	aview1 offset: 0
       	aview1 stride: 2
       	aview1 data:	0   2   4   6   8   

       	aview2 data has odd values:
       	aview2 num elements: 5
       	aview2 offset: 1
       	aview2 stride: 2
       	aview2 data:	1   3   5   7   9   

  After destroyViewAndData(aview1) call
       	Num views in group A: 1
       	Num buffers in datastore: 1
       	Buffer num elements: 10
       	aview2 data still has its odd values:	1   3   5   7   9 

-----------------------------------------------
Example 3: Two Views and One Buffer (View Copy)
-----------------------------------------------

The third example illustrates another  Sidre usage pattern in which multiple 
views share the same data in a single buffer. Before we start, we verify that 
the datastore contains one buffer and that the "A" group contains one view. 
This is the Sidre buffer, group, and view state that exists 
at the end of the second example.

We begin by making a copy of the view in group "A" in group "B". The new 
view is identical in description and data associated with it as the original
view. In particular, the data in the new view *is the same data* in the
original view. Recall that Sidre copy operations for groups and views are
**shallow copy operations**. We verify this, by printing the values of the
array associated with each view and also the base address of each array.

Next, we destroy the "A" group which owned the original view. We verify that
the view copy remains in the "B" group and its data is still intact. When
we destroyed the "A" group, its view is *detached* but not destroyed. We see
that it is still exist and is attached to the buffer.

Finally, we destroy the "B" group. Similar to the destruction of the "A" group,
the view is detached from the "B" group. In particular, the view remains 
attached to the buffer, which still exists in the datastore and remains 
allocated.

The last two operations in this example are intended to emphasize the 
explicit nature of Sidre operations. For example, when a group is destroyed,
its views and their data remain intact. There are specific methods provided
to destroy views and to deallocate and destroy buffers. This was an early
design choice for Sidre to yield maximum flexibility in defining data 
hierarchies while keeping internal bookkeeping simple. The example following 
this one continues on this point.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex3_twoviews_onebuffer_copy_start 
   :end-before: _ex3_twoviews_onebuffer_copy_end
   :language: C++

The output of the code is::

  Datastore start state
       	Num buffers in datastore: 1
       	Num views in group A: 1
       	Num views in group B: 0

  After copying aview2 to group B
       	Num buffers in datastore: 1
       	Num views in group A: 1
       	Num views in group B: 1
       	aview2 in A group has values:	1   3   5   7   9   
       	aview2 in B group has values:	1   3   5   7   9   

       	Base address of array in A group: 0xc595b8
       	Base address of array in A group: 0xc595b8

  After destroyGroup(A_grp) call:
       	Num buffers in datastore: 1
        Is buffer allocated? 1 
        Num views attached to buffer: 2
       	Num views in group B: 1
       	aview2 in B group has values:	1   3   5   7   9

  After destroyGroup(B_grp) call:
        Num buffers in datastore: 1
        Is buffer allocated? 1
        Num views attached to buffer: 2
  
-------------------------------------
Example 4: More Basic Mechanics
-------------------------------------

The last example should help to make clear the point made at the beginning of
this section end the end of example three about all Sidre operations being 
simple and explicit.

We start with a fresh datastore with one group "A" in the root group. We
create a view "aview" in the group and allocate it. The output of the 
code below shows that everything is working as expected. We know what to expect
because this part is the same as the first part of example one.

We destroy the view and we see that the view is gone from the group, but
that the buffer is still in the datastore and allocated. This is so because
we did not deallocate the view or destroy the buffer.

Next, we recreate the view, attach the buffer to it and describe it as before. 
This restores the state of everything before we destroyed the view.

Then, we deallocate the buffer and verify that the buffer and view are no 
longer allocated, but are still described. Also, the view is still attached 
to the buffer.

Lastly, we destroy the buffer. It has been deallocated and gone from the 
datastore.  

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex4_more_mechanics_start
   :end-before: _ex4_more_mechanics_end
   :language: C++

The output of the code is::

  After A_grp->createViewAndAllocate() call:
        Num views in group A: 1
        Num elements in view: 10
        Is view allocated? 1
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After A_grp->destroyView() call:
        Num views in group A: 0
        Num buffers in datastore: 1
        Num views attached to buffer: 0
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After recreating view, attaching buffer, and describing data:
        Num views in group A: 1
        Num elements in view: 10
        Is view allocated? 1
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 1

  After buffer deallocate call:
        Num elements in view: 10
        Is view allocated? 0
        Num buffers in datastore: 1
        Num views attached to buffer: 1
        Num elements in buffer array: 10
        Is buffer allocated? 0

  After destroy buffer call:
        Num elements in view: 10
        Is view allocated? 0
        Num buffers in datastore: 0

The point of this final example is to demonstrate that most Sidre operations are
atomic and do precisely what you expect with few side effects. If this were not
the case, we believe that the Sidre API would become bloated and confusing to 
be able to support the range of usage patterns applied by Axom users.

