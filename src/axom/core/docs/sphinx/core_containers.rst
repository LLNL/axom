.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core containers
******************************************************

.. warning::
   Axom's containers are currently the target of a refactoring effort. The 
   following information is subject to change, as are the interfaces for the 
   containers themselves. When changes are made to the Axom containers, the
   changes will be reflected here.

Axom Core contains the ``Array``, ``ArrayView``, and ``StackArray`` classes.  
Among other things, these data containers facilitate porting code that uses 
``std::vector`` to GPUs.

``Array`` is a multidimensional contiguous container template. In the 
1-dimensional case, this class behaves similar to ``std::vector``. In higher 
dimensions, some vector-like functionality, such as ``push_back``, are not
available and multidimensional-specific operations will mirror the 
``ndarray`` provided by ``numpy`` when possible. 

The ``Array`` object manages all memory. Typically, the ``Array`` object will
allocate extra space to facilitate the insertion of new elements and minimize the number of reallocations.
The actual capacity of the array (i.e., total number of elements that
the ``Array`` can hold) can be queried via the ``capacity()`` method.
When allocated memory is used up, inserting a new element triggers a
reallocation.  At each reallocation, extra space is allocated
according to the **resize_ratio** parameter, which is set to 2.0
by default. To return all extra memory, an application can call
`shrink()`.

.. warning:: Reallocations tend to be costly operations in terms of performance.
  Use ``reserve()`` when the number of nodes is known a priori, or
  use a constructor that takes an actual size and capacity when possible.
  
.. note:: The Array destructor deallocates and returns all memory associated
  with it to the system.

Here's an example showing how to use ``Array`` instead of ``std::vector``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraybasic_start
   :end-before: _arraybasic_end
   :language: C++

The output of this example is::

  Length of a = 3
  After appending a value, a's length = 4
  Array a = [2, 5, 11, 4]
  After inserting two values, Array a = [2, 5, 6, 11, 1, 4]

In the future, it will be 
possible to take "views" of the underlying array data that allow for 
flexible reinterpretation via different striding.

Applications commonly store tuples of data in a flat array or a ``std::vector``.
In this sense, it can be thought of as a two-dimensional array.  ``Array``
supports arbitrary dimensionalities but an alias, ``MCArray``
(or Multi-Component Array), is provided for this two-dimensional case.

The ``MCArray`` (i.e., ``Array<T, 2>``) class formalizes tuple storage, as shown in the next example.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraytuple_start
   :end-before: _arraytuple_end
   :language: C++

The output of this example is::

  MCArray b with 2 3-tuples = [
    [1, 4, 2]
    [8, 0, -1]
  ]
  MCArray b with 4 3-tuples = [
    [1, 4, 2]
    [0, -1, 1]
    [1, -1, 0]
    [8, 0, -1]
  ]

It is also often useful to wrap an external, user-supplied buffer without taking ownership
of the data.  For this purpose Axom provides the ``ArrayView`` class, which is a lightweight
wrapper over a buffer that provides one- or multi-dimensional indexing/reshaping semantics.
For example, it might be useful to reinterpret a flat (one-dimensional) array as a two-dimensional array.
This is accomplished via ``MCArrayView`` which, similar to the ``MCArray`` alias, is an alias for ``ArrayView<T, 2>``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _extbuffer_start
   :end-before: _extbuffer_end
   :language: C++

The output of this example is::

  Array a = [2, 5, 6, 11, 1, 4]
  MCArrayView c with 3 2-tuples = [
    [2, 5]
    [6, 11]
    [1, 4]
  ]
  Array a and MCArrayView c use the same memory, a's internal buffer.
  Array a = [1, 5, 6, 9, 1, 4]
  MCArrayView c with 3 2-tuples = [
    [1, 5]
    [6, 9]
    [1, 4]
  ]

.. note:: The set of permissible operations on an ``ArrayView`` is somewhat limited,
  as operations that would cause the buffer to resize are not permitted.

The ``StackArray`` class is a work-around for a limitation in older versions
of the nvcc compiler, which do not capture arrays on the stack in device 
lambdas.  More details are in the API documentation and in the tests.
