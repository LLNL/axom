.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core acceleration
******************************************************

Axom Core provides support to ease the implementation of user applications that
intend to support execution utilizing accelerators. In order to supply this,
Axom provides a simple API to access memory and computation, using RAJA
and Umpire to access accelerators.

.. warning:: 
   You do need to compile with RAJA and Umpire to fully leverage this Axom functionality. 
   While both submodules will operate without these tools, they will offer little.

The memory management submodule offers a singular interface which allows the user to leverage either C++ memory functions or Umpire, 
depending on the availability of Umpire at compilation. It supports a simple operation set: allocate, deallocate, reallocate, and copy. 
If Umpire is in use, an allocator can be specified -- see Umpire documentation(make this a link) for more details. Note that if an 
Umpire-reliant function is used when Umpire is not present, such as for allocator setting, the function fails without a fault, so one 
piece of code can handle standard C++ or C++ with Umpire. 

Here is an example of using Axom’s memory management tools:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _membasic_start
   :end-before: _membasic_end
   :language: C++

Throughout Axom, acceleration is increasingly supported. Both internally, and to support users, Axom Core offers an 
interface that, using RAJA and Umpire internally, provides easy access to for-loop level acceleration via the for-all model, 
which applies a given lambda function for every index in range.  

Here is an example of Axom in motion with basic sequential execution:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _exebasic_start
   :end-before: _exebasic_end
   :language: C++

Here's the same loop from the above snippet, this time with CUDA:

.. literalinclude:: ../../examples/core_acceleration.cpp
   :start-after: _cudaexebasic_start
   :end-before: _cudaexebasic_end
   :language: C++

For more advanced functionality, see the RAJA documentation(link here), and the Umpire documentation(link here). Axom 
allows for these to be used in lieu of Axom’s API.  
