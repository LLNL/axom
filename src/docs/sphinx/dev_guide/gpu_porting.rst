.. ## Copyright (c) 2017-2022, Lawrence Livermore National Security,con LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _gpu-porting-label:

*******************************
GPU Porting in Axom
*******************************

Axom uses RAJA and Umpire as the main workhorses for GPU porting,
with a set of convenience macros and `axom`-namespaced wrappers encapsulating
commonly used RAJA and Umpire functions, and preset execution spaces for
host/device execution.

For the user's guide on using GPU utilities, see also
:ref:`Core Acceleration<core-acceleration>`


.. _gpu-macros-label:

===============
Macros 
===============

Axom's macros can be found in file
`axom/core/Macros.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/Macros.hpp>`_

Most of the GPU-related macros are used to guard device code for compilation
time or for ``__host__ __device__`` decoration of functions/lambdas.

For guarding device code:

.. literalinclude:: ../../../axom/core/Macros.hpp
   :start-after:  _guarding_macros_start
   :end-before:  _guarding_macros_end
   :language: C++

For ``__host__`` or `` __device__`` decoration, or both:

.. literalinclude:: ../../../axom/core/Macros.hpp
   :start-after:  _decorating_macros_start
   :end-before:  _decorating_macros_end
   :language: C++

.. _gpu-memory-label:

===============
Memory 
===============

Axom's memory management routines can be found in the file
`axom/core/memory_management.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/memory_management.hpp>`_


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Memory Management Routines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Umpire has the concept of "allocators" associated with each
`memory resource type`_ (``umpire::resource::MemoryResourceType``)

To allocate memory on a particular resource, you use the ID for the allocator
associated with the MemoryResourceType.

You are able to set a default allocator, whereby all your memory allocations
will go on the resource associated with the allocator unless otherwise
specified::

  //---------------------------------------------------------------------------
  // Getters and Setters for Allocators and ID
  //---------------------------------------------------------------------------

  // Gets the allocator ID associated with the given MemoryResourceType
  int axom::getUmpireResourceAllocatorID(umpire::resource::MemoryResourceType resource_type)

  // Sets the default allocator associated with the given MemoryResourceType
  void axom::setDefaultAllocator(umpire::resource::MemoryResourceType resource_type)

  // Sets the default allocator associated with the given allocator ID
  void axom::setDefaultAllocator(int allocId)

  // Returns the current default allocator ID
  int axom::getDefaultAllocatorID()

  //---------------------------------------------------------------------------
  // Memory Allocation Functions
  //---------------------------------------------------------------------------

  // Allocates n number of T elements in the memory space with the given
  // allocator ID. If no ID provided, uses the current default
  T* axom::allocate(std::size_t n, int allocID = getDefaultAllocatorID())

  //  Deallocates the given pointer
  void axom::deallocate(T*& p)

  // Reallocates pointer to given size and with the given allocator ID.
  // If no ID provided, uses the current default.
  T* axom::reallocate(T* p, std::size_t n, int allocID = getDefaultAllocatorID())

  // Copies numbytes from src to dst
  axom::copy(void* dst, const void* src, std::size_t numbytes)

.. note::

  When Axom is built without Umpire, the getters and setters shown above
  become no-ops or are undefined, while the memory allocation functions
  default to C++ standard library functions with only allocation on the
  host (CPU):

  * ``axom::allocate`` calls ``std::malloc``
  * ``axom::deallocate`` calls ``std::free``
  * ``axom::reallocate`` calls ``std::realloc``
  * ``axom::copy`` calls ``std::memcpy``


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
MemorySpace
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. literalinclude:: ../../../axom/core/memory_management.hpp
   :start-after:  _memory_space_start
   :end-before:  _memory_space_end
   :language: C++

These ``axom::MemorySpace`` enum values are used to define where in memory
``axom::Array`` and ``axom::ArrayView`` types are defined.

``Dynamic`` allows you to define the location at runtime, with some caveats
(see :ref:`Core Containers<core-containers>` for more details and examples).


^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Useful Links
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`Umpire Tutorial`_ - 1st two sections cover Allocators and Resources


.. _gpu-kernels-label:

===============
Kernels
===============

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
axom::for_all
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``axom::for_all`` can be found in
`axom/core/execution/for_all.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/execution/for_all.hpp>`_

``axom::for_all`` is a wrapper around `RAJA forall`_, a simple for loop.

This is used the most in Axom, from unit tests to example drivers for loops
on device::

  template <typename ExecSpace, typename KernelType>
  void axom::for_all(const IndexType& N, KernelType&& kernel)

  template <typename ExecSpace, typename KernelType>
  void axom::for_all(const IndexType& begin, const IndexType& end, KernelType&& kernel)

.. note::

  When Axom is built without RAJA, ``axom::for_all`` becomes a for loop on
  host (CPU).

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
RAJA::kernel
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
``RAJA::kernel`` is for more complex and/or nested loops.

This is used infrequently, mainly seen only in a `few unit tests`_.

Your general go-to will be ``axom::for_all``.

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Useful Links
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
`RAJA Loops`_ - Covers for_all, kernel, and other loops


.. _gpu-execution-spaces-label:

================================
Execution Spaces & Policies
================================

Axom's execution spaces can be found in
`axom/core/execution/execution_space.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/core/execution/execution_space.hpp>`_

Axom's execution spaces are derived from an execution_space traits class,
containing RAJA execution policies and default Umpire memory allocators
associated with each space.

Axom currently supports `4 execution spaces`_:

* ``SEQ_EXEC`` - Sequential execution policies on host
* ``OMP_EXEC`` -  OpenMP execution policies on host
* ``CUDA_EXEC`` - CUDA execution policies in Unified Memory (host + device)
* ``HIP_EXEC`` - HIP execution policies in Unified Memory (host + device)

Additionally, ``HIP_EXEC`` and ``CUDA_EXEC`` are templated by the
number of threads and SYNCHRONOUS or ASYNC execution.

.. literalinclude:: ../../../axom/core/execution/internal/cuda_exec.hpp
   :start-after:  _cuda_exec_start
   :end-before:  _cuda_exec_end
   :language: C++

Each execution space (``axom::execution_space<ExecSpace>``) provides:

* Axom policies that are type aliases of RAJA policies to be used with kernels,
  RAJA types, and RAJA operations

  * ``loop_policy`` - For `RAJA scans`_ and other operations; ``axom::for_all``
    uses the loop_policy from the templated execution space.
  * ``reduce_policy`` - For  `RAJA reduction types`_
  * ``atomic_policy`` - For `RAJA atomic operations`_
  * ``sync_policy`` - For Axom's `synchronize`_ function, which is a wrapper
    around ``RAJA::synchronize()``

* Umpire allocator defaults

  * ``memory_space`` - The memory space abstraction for use by
    :ref:`Core Containers<core-containers>` like ``axom::Array``.
  * ``allocatorID()`` - Gets the allocator ID for the Umpire resource to use in
    this execution space.

* General information on the execution space

  * ``name()`` - Name of the execution space
  * ``onDevice()`` - Is the execution space on device? (True/False)
  * ``valid()`` - Is the execution space valid? (True)
  * ``async()`` - Is the execution space asynchronous? (True/False)

The mint component also provides a set of nested execution policies
located at
`axom/mint/execution/internal/structured_exec.hpp <https://github.com/LLNL/axom/blob/develop/src/axom/mint/execution/internal/structured_exec.hpp>`_
to be used with
``RAJA::kernel`` e.g. for iterating over mint meshes.

.. note::

  When Axom is built without RAJA, only ``SEQ_EXEC`` is available
  for host (CPU) execution. When Axom is built without Umpire, only
  ``SEQ_EXEC`` and ``OMP_EXEC`` is available for host (CPU) execution.


.. _gpu-tips-label:

================================
General, Rough Porting Tips
================================

* Start with figuring out what memory you need on device, and use Axom's Array
  class and memory_managment functions to do the allocations.
* Using a barebones Axom's for_all kernel with a device policy, attempt to
  access/ manipulate the memory on device.

  * This does not have to be logically correct, you are just making sure you
    have done the allocation correctly (e.g. not segfaulting).

* Add the functions you want to call on device to the ``axom::for_all`` kernel.

  * Does not have to be logically correct.

* When you attempt to recompile your code, unless you have already done the
  ``__host__ __device__`` decoration correctly already, the compiler should
  complain about what functions are not yet compiling on device, and that you
  need to add ``__host__ __device__`` decorators::

    error: calling a __host__ function("XXX") from a __host__ __device__ function("YYY") is not allowed

  * Even after you add the decorators to the functions the compiler complains
    about, it is likely recompiling will introduce complaints about more
    functions (the functions being the non-decorated functions your
    newly-decorated functions are calling). You will want to keep decorating
    until all the complaints are gone.
  * Most of the C++ standard library is not available on device. Your options
    are Axom's equivalent functions/classes if it exists, or to add your
    own / rewrite the code to not use standard library.

* With no more decorating complaints from the compiler, write the logically
  correct kernel.

  * If at this point your kernel is not working/segfaulting, it's hopefully a
    logical error, and you can debug the kernel without diving into debugging
    tools.
  * Utilize ``printf()`` for debugging output
  * Try using the ``SEQ_EXEC`` execution space

.. _memory resource type: https://github.com/LLNL/Umpire/blob/develop/src/umpire/resource/MemoryResourceTypes.hpp#L63
.. _Umpire Tutorial: https://umpire.readthedocs.io/en/develop/sphinx/tutorial.html
.. _RAJA forall: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/loop_basic.html#simple-loops-raja-forall
.. _few unit tests: https://github.com/LLNL/axom/search?q=RAJA%3A%3Akernel
.. _RAJA Loops: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/loop_basic.html#elements-of-loop-execution
.. _4 execution spaces: https://github.com/LLNL/axom/tree/develop/src/axom/core/execution/internal
.. _RAJA scans: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/scan.html#scans
.. _RAJA reduction types: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/reduction.html#reduction-operations
.. _RAJA atomic operations: https://raja.readthedocs.io/en/develop/sphinx/user_guide/feature/atomic.html#atomics
.. _synchronize: https://github.com/LLNL/axom/blob/482749b322057ca7ef9a6d786e948692d9f72246/src/axom/core/execution/synchronize.hpp
