.. ## Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Core utilities
******************************************************

The `axom::utilities` namespace contains basic useful functions.  Often
these have started out in another, higher-level library and proven so
useful they're "promoted" to Axom Core.  In some cases, `axom::utilities`
brings functionality from recent standards of C++ to platforms restricted
to older compilers.

Here a function showing the usage of string and filesystem facilities available
in Axom Core.

.. literalinclude:: ../../examples/core_utilities.cpp
   :start-after: _fs_string_start
   :end-before: _fs_string_end
   :language: C++

Axom Core also includes a `Timer` class.  Here, we time the preceding
filesystem example snippet.

.. literalinclude:: ../../examples/core_utilities.cpp
   :start-after: _timer_start
   :end-before: _timer_end
   :language: C++

