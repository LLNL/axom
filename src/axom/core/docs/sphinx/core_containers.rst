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
1-dimensional case, this class behave similar to ``std::vector``. In higher 
dimensions, some vector-like functionality will be unavailable, e.g., 
``push_back``, and multidimensional-specific operations will mirror the 
``ndarray`` provided by ``numpy`` when possible. In the future, it will be 
possible to take "views" of the underlying array data that allow for 
flexible reinterpretation via different striding.

``MCArray`` (or Multi-Component Array) is a container template for a contiguous
array of tuples. In this sense, it can be thought of as a two-dimensional array.``MCArray`` will be removed in the future as ``Array`` will provide all the 
this multidimensional functionality of ``MCArray``..

Here's an example showing how to use ``MCArray`` instead of ``std::vector``.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraybasic_start
   :end-before: _arraybasic_end
   :language: C++

Applications commonly store tuples of data in a flat array or a ``std::vector``.
The ``MCArray`` class formalizes tuple storage, as shown in the next example.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraytuple_start
   :end-before: _arraytuple_end
   :language: C++

The ``MCArray`` class can use an external memory buffer, with the restriction 
that operations that would cause the buffer to resize are not permitted.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _extbuffer_start
   :end-before: _extbuffer_end
   :language: C++

Similar to ``std::vector``, when using an internal array the ``Array`` and 
``MCArray`` classes can reserve memory in anticipation of future growth as 
well as shrink to the memory in use.

The ``StackArray`` class is a work-around for a limitation in older versions
of the nvcc compiler, which do not capture arrays on the stack in device 
lambdas.  More details are in the API documentation and in the tests.
