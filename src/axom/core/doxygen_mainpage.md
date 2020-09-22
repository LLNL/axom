Core {#coretop}
=============================

Shared utility functions and classes supporting all Axom components.
These hide the difference between versions of C++ and between supported computing platforms.

- [Utilities](@ref axom::utilities) include
  - A [Timer](@ref axom::utilities::Timer) class to measure elapsed time
  - Templated arithmetic and logical comparison functions
  - A processAbort() function
- [Numerics](@ref axom::numerics) include
  - A [Matrix](@ref axom::numerics::Matrix) class and operators
  - Polynomial and linear system solvers
