.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core containers
******************************************************

Axom Core contains the Array and StackArray classes.  Among other things, these
data containers facilitate porting code that uses `std::vector` to the GPU.


Here's an example showing how to use Array instead of `std::vector`.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraybasic_start
   :end-before: _arraybasic_end
   :language: C++

Applications commonly store tuples of data in a flat array or a `std::vector`.
The Array class formalizes tuple storage, as shown in the next example.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _arraytuple_start
   :end-before: _arraytuple_end
   :language: C++

The Array class can use an external memory buffer, with the restriction that 
operations that would cause the buffer to resize are not permitted.

.. literalinclude:: ../../examples/core_containers.cpp
   :start-after: _extbuffer_start
   :end-before: _extbuffer_end
   :language: C++

Similar to `std::vector`, when using an internal array the Array class can
reserve memory in anticipation of future growth as well as shrink to just the
memory currently in use.

The StackArray class is a work-around for a limitation in the nvcc compiler,
which can't capture arrays on the stack in device lambdas.  More details are in
the API documentation and in the tests.
