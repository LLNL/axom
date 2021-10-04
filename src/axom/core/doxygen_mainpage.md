Core {#coretop}
=============================

Shared utility functions and classes supporting all Axom components.
These hide the difference between versions of C++ and between supported computing platforms.

- [Utilities](@ref axom::utilities) include
  - A [Timer](@ref axom::utilities::Timer) class to measure elapsed time
  - [File system operations](@ref axom::utilities::filesystem)
  - [String manipulation](@ref axom::utilities::string) 
  - Routines for getting information about the machine and user running the code (axom::utilities::getHostName() and axom::utilities::getUserName() methods) and gracefully exiting execution (processAbout() method)
  - Routines for getting information about the Axom version (axom::getVersion() method) and how Axom was configured and built (axom::about() method)
  - Templated arithmetic and logical comparison functions
- [Containers]
  - [Basic linear array](@ref axom::Array)
  - [Multi-component array](@ref axom::MCArray)
  - [Stack array](@ref axom::StackArray) for passing stack arrays to device kernels (e.g., GPU) via lambda capture properly.
- [Numerics](@ref axom::numerics) include
  - A [Matrix](@ref axom::numerics::Matrix) class and operators
  - Polynomial and linear system solvers
