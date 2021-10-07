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
and how the description of data in a view is not completely dependent on the
data in the corresponding buffer. The code snippets shown and described here
exist in the file ``axom/src/axom/sidre/examples/sidre_data_vs_metadata.cpp``,
which can be built and run to experiment with if you wish.

The starting point for the examples is a simple Sidre datastore group hierarchy
in which the root group has two child groups named "A" and "B". This structure
is generated as follows:

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _datastore_initial_start
   :end-before: _datastore_initial_end
   :language: C++

At this point, the datastore has zero buffers and the groups have no views.

--------------------------------
Example 1: One View, One Buffer 
--------------------------------

In the first exaample, we first create a view in group "A" describing an integer
array of length 10 and allocate it. This creates a buffer in the datastore
which owns the array. Then, we get a handle to the buffer from the view
from which we get a pointer to the start of the array and initialize the 
array values. To give some insight into the internal Sidre mechanics, we access
and print various pieces of information showing the number of buffers in the 
datastore, the number of views in group "A", the number of elements in the view,
the number of views attached to the buffer, the number of elements in the 
the buffer, and the buffer array value at slot 5.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_create_start
   :end-before: _ex1_oneview_onebuffer_create_end
   :language: C++

The output printed by the code is::

  After A_grp->createViewAndAllocate() call
        Num buffers in datastore: 1
        Num views in group A: 1
        Num elements in view: 10
        Num views attached to buffer: 1
        Num elements in buffer array: 10

  After initialization of view array
        Value of elt 5 in buffer array (expect 7): 7

Next, we deallocate the view and show that its decription remains intact; 
for example, num elements is 10. But, the view itself is no longer allocated
meaning that the data is no longer accessible via the view. The allocated 
buffer still exists in the datastore. Then, we allocate the buffer (again). 
Since we have not changed the view description and the view is still attached
to the buffer (we did not explicitly detach it), the allocate operation 
re-associates the buffer data with the view and does not create a new buffer 
in the datastore. We verify this by showing that the array value at slot 5 is 
the same as before and that the number of buffers in the datastore is still one.

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex1_oneview_onebuffer_deallocalloc_start
   :end-before: _ex1_oneview_onebuffer_deallocalloc_end
   :language: C++

The output of the code is::

  After view deallocate call, the data is no longer accessible
  via the view, but its description remains.
        Is view allocated? 0
        Num elements in view: 10
        Num buffers in datastore: 1
        Is buffer allocated? 0

  After allocating view again...
        Is view allocated? 1
        Num elements in view: 10 
        Value of elt 5 in view array (expect 7): 7
        Num buffers in datastore: 1

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

----------------------------
Example 2
----------------------------

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex2_twoviews_onebuffer_start
   :end-before: _ex2_twoviews_onebuffer_end
   :language: C++

----------------------------
Example 3
----------------------------

.. literalinclude:: ../../examples/sidre_data_vs_metadata.cpp
   :start-after: _ex3_twoviews_onebuffer_copy_start 
   :end-before: _ex3_twoviews_onebuffer_copy_end
   :language: C++

