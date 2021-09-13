.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core acceleration
******************************************************

Axom's `core` component provides several utilities for user applications that intend
to support execution on hardware accelerators.
Axom lets users control execution space using RAJA and memory space using Umpire.
As noted in the `RAJA documentation <https://raja.readthedocs.io/en/main/index.html>`_,
developers of high-performance computing applications have many options for running
code: on the CPU, using OpenMP, or on GPU hardware accelerators, and the options are
constantly developing.  RAJA controls *where* code runs; Umpire moves data between
memory spaces.

.. note::
   Axom's memory management and execution space APIs have default implementations when Axom is not
   configured with Umpire and RAJA, respectively, and can be used in all Axom configurations.

The memory management API allows the user to leverage either C++ memory functions or Umpire, 
depending on the availability of Umpire at compilation. It supports a simple operation set: allocate, deallocate, reallocate, and copy. 
If Umpire is in use, an allocator can be specified -- see `Umpire documentation <https://readthedocs.org/projects/umpire>`_ 
for more details. Note that the fall-back to C++ memory functions is automatic, so the same piece of code can handle standard C++ 
or C++ with Umpire. However, to use advanced features, such as accessing unified memory, Umpire must be enabled, otherwise errors 
will occur at compilation.

Here is an example of using Axom's memory management tools:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _membasic_start
   :end-before: _membasic_end
   :language: C++

Throughout Axom, acceleration is increasingly supported. Both internally, and to support users, Axom Core offers an 
interface that, using RAJA and Umpire internally, provides easy access to for-loop level acceleration via the parallel-for idiom, 
which applies a given lambda function for every index in range.  

Here is an Axom example showing sequential execution:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _exebasic_start
   :end-before: _exebasic_end
   :language: C++

Here's the same loop from the above snippet, this time with CUDA:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _cudaexebasic_start
   :end-before: _cudaexebasic_end
   :language: C++

For more advanced functionality, users can directly call RAJA and Umpire.
See the `RAJA documentation <https://raja.readthedocs.io>`_
and the `Umpire documentation <https://umpire.readthedocs.io>`_.
