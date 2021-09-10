.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core acceleration
******************************************************

Axom Core provides support for user applications that target accelerators.
Axom lets users control execution space using RAJA and memory space using Umpire.

.. warning:: 
   You do need to compile with RAJA and Umpire to fully leverage this Axom functionality. 
   While both memory management and execution spaces in Axom will operate without these 
   tools, they will offer little.

The memory management facility offers a singular interface which allows the user to leverage either C++ memory functions or Umpire, 
depending on the availability of Umpire at compilation. It supports a simple operation set: allocate, deallocate, reallocate, and copy. 
If Umpire is in use, an allocator can be specified -- see `Umpire documentation <https://readthedocs.org/projects/umpire>`_ 
for more details. Note that the fall-back to C++ memory functions is automatic, so the same piece of code can handle standard C++ 
or C++ with Umpire. However, to use advanced features, such as accessing unified memory, Umpire must be enabled, otherwise errors 
will occur at compilation.

Here is an example of using Axom’s memory management tools:

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

For more advanced functionality, see the `RAJA documentation <https://raja.readthedocs.io>`_, and the `Umpire documentation <https://umpire.readthedocs.io>`_. Axom 
allows for these to be used in lieu of Axom’s API.  
