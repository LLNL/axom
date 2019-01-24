.. ##
.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

***********************************
Surface mesh point queries: C++ API
***********************************

Codes written in C++ may use the object-oriented C++ APIs to perform in/out and
signed distance queries.  In addition to language choice, the C++ API lets a code work
with more than one mesh at the same time.  Unlike the C API, the C++ API for
in/out and signed distance queries has no initializer taking a file name:
readying the mesh is a separate, prior step.

In/out Octree
-------------

The C++ in/out query is provided by the ``quest::InOutOctree`` class, from the
following header.  See axom/src/axom/quest/tests/quest_inout_octree.cpp.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_include_start
   :end-before: _quest_inout_cpp_include_end
   :language: C++

Some typedefs are useful for the sake of brevity.  The class is templated on
the dimensionality of the mesh.  Currently, only 3D meshes are supported;
here ``DIM`` equals 3.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_typedef_start
   :end-before: _quest_inout_cpp_typedef_end
   :language: C++

Instantiate the object using ``GeometricBoundingBox bbox`` and a mesh, and
generate the index.

.. literalinclude:: ../../tests/quest_inout_octree.cpp
   :start-after: _quest_inout_cpp_init_start
   :end-before: _quest_inout_cpp_init_end
   :language: C++

Test a query point.
::

   SpacePt pt = SpacePt::make_point(2., 3., 1.);
   bool inside = octree.within(pt);

All cleanup happens when the index object's destructor is called 
(in this case, because the variable goes out of scope).

Octree spatial data index
=========================

The ``quest::InOutOctree`` class implements the in/out point query against a surface mesh.
This class is a subclass of the ``quest::SpatialOctree``, which can be used by client codes
to implement custom spatial query acceleration data structures.  An N-dimensional
octree spatial index divides its region of interest into blocks by intersecting a
subdivision point with N perpendicular planes, resulting in 2^N child blocks contained
in the parent block.  A common usage is in 3D, hence the name octree.

The ``quest::SpatialOctree`` is templated on dimensionality and block type.  The block type
can contain any metadata needed to characterize the spatial extent it represents.  The
reference block type ``quest::BlockData``, for example, stores a numeric ID (possibly
the ID of the partition point).  The ``InOutOctree`` uses a class for its block type 
that stores data for the in/out query along with the partition point ID.  To be used 
as a block data type, a class must be default-constructable and provide the functions 
``bool isLeaf()``, ``void setInternal()``, ``void setNonBlock()``, and ``bool isBlock()``.

Signed Distance
---------------

The C++ signed distance query is provided by the ``quest::SignedDistance`` class,
which wraps an instance of ``primal::BVHTree``.
Examples from axom/src/axom/quest/tests/quest_signed_distance.cpp.
Class header:

.. literalinclude:: ../../tests/quest_signed_distance.cpp
   :start-after: _quest_distance_cpp_include_start
   :end-before: _quest_distance_cpp_include_end
   :language: C++

The constructor takes several arguments.  Here, ``surface_mesh`` is a pointer to
a triangle surface mesh.  The second argument indicates the mesh is watertight.
The third and fourth arguments are used to build the underlying BVH tree
spatial index and indicate that BVH tree buckets will be split after 25 objects
and that the BVH tree will contain at most 10 levels.  These are safe default 
values and can be adjusted after benchmarking shows the need.  Note that the
second and subsequent arguments to the constructor correspond to the 
``quest::signed_distance_set`` functions in the C API.

As with the ``InOutOctree``, the class is templated on the dimensionality
of the mesh, with only 3D meshes being supported.

.. literalinclude:: ../../tests/quest_signed_distance.cpp
   :start-after: _quest_distance_cpp_init_start
   :end-before: _quest_distance_cpp_init_end
   :language: C++

Test a query point.
::

   axom::primal::Point< double,3 > pt =
     axom::primal::Point< double,3 >::make_point(2., 3., 1.);
   double signedDistance = signed_distance.computeDistance(pt);

The object destructor takes care of all cleanup.
