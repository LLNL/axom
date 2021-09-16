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

Axom Core contains the ``Array``, ``MCArray``, and ``StackArray`` classes.  
Among other things, these data containers facilitate porting code that uses 
``std::vector`` to GPUs.

``Array`` is a multidimensional contiguous container template. In the 
1-dimensional case, this class behaves similar to ``std::vector``. In higher 
dimensions, some vector-like functionality, such as ``push_back``, are not
available and multidimensional-specific operations will mirror the 
``ndarray`` provided by ``numpy`` when possible. In the future, it will be 
possible to take "views" of the underlying array data that allow for 
flexible reinterpretation via different striding.

``MCArray`` (or Multi-Component Array) is a container template for a contiguous
array of tuples. In this sense, it can be thought of as a two-dimensional array.
``MCArray`` will be removed in the future as ``Array`` will provide all the 
multidimensional functionality of ``MCArray``..

Here's an example showing how to use ``MCArray`` instead of ``std::vector``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraybasic_start
   :end-before: _arraybasic_end
   :language: C++

The output of this example is::

  Length of a = 3
  After appending a value, a's length = 4
  MCArray a = [2, 5, 11, 4]
  After inserting two values, MCArray a = [2, 5, 6, 11, 1, 4]

Applications commonly store tuples of data in a flat array or a ``std::vector``.
The ``MCArray`` class formalizes tuple storage, as shown in the next example.

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

The ``MCArray`` class can use an external memory buffer, with the restriction 
that operations that would cause the buffer to resize are not permitted.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _extbuffer_start
   :end-before: _extbuffer_end
   :language: C++

The output of this example is::

  MCArray a = [2, 5, 6, 11, 1, 4]
  MCArray c with 3 2-tuples = [
    [2, 5]
    [6, 11]
    [1, 4]
  ]
  MCArrays a and c use the same memory, a's internal buffer.
  MCArray a = [1, 5, 6, 9, 1, 4]
  MCArray c with 3 2-tuples = [
    [1, 5]
    [6, 9]
    [1, 4]
  ]

Similar to ``std::vector``, when using an internal array the ``Array`` and 
``MCArray`` classes can reserve memory in anticipation of future growth as
well as shrink to the memory in use. Please refer to the 
`Array API Documentation <../../../../doxygen/html/classaxom_1_1Array.html>`_ 
and 
`MCArray API Documentation <../../../../doxygen/html/classaxom_1_1Array.html>`_ 
for details.

The ``StackArray`` class is a work-around for a limitation in older versions
of the nvcc compiler, which do not capture arrays on the stack in device 
lambdas.  More details are in the API documentation and in the tests.
