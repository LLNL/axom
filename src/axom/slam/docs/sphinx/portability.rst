.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _portability-label:

Host and device use
===================

Slam supports host construction and device access for selected configurations.
Check the complete type and the operations a kernel calls, not just its
indirection policy. A map backed by ``ArrayView`` still refers to a set, and a
submap still refers to its parent map.

Before using an object in a kernel, check that:

* Every operation on the kernel's call path is device-callable with the chosen
  compiler and standard library.
* Referenced sets, parent objects and buffers reside in memory the device can access.
* Those objects and allocations remain valid until the kernel finishes,
  including asynchronous execution.

``TriviallyCopyableRepresentation<T>`` checks only the C++ object
representation. The virtual interfaces are not a portable cross-device dispatch
mechanism.

C++20 ranges integration
------------------------

Host code that relies on Slam's standard-ranges customizations must include the
host-only header explicitly:

.. code-block:: cpp

   #include "axom/slam/Ranges.hpp"

This header marks ``GenericRangeSet`` configurations with ``NoIndirection``
and ``NoSubset`` as borrowed ranges. This includes ``RangeSet`` and
``PositionSet``. Their iterators store the range state by value and can outlive
the range object. This does not make arbitrary indirection sets, maps or
submaps borrowed ranges. Their iterators can depend on other objects or buffers.

``Ranges.hpp`` is installed with Slam but excluded from the unified
``axom/slam.hpp`` header. Device-facing translation units therefore do not
acquire a dependency on ``<ranges>`` through that header. Slam's standard-ranges
integration is for host code.

Optional results in kernels
---------------------------

Some Slam APIs return ``std::optional``. Device use depends on the compiler,
standard library and build flags, so compile and exercise the actual accessor
path in the target configuration.

On a supported device path, check ``has_value()`` before dereferencing with
``operator*``, or use ``value_or`` with a device-compatible value type. Avoid
``value()`` in kernels because an empty optional can invoke a throwing,
host-only helper.
