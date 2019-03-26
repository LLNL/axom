# Axom Software Release Notes

Notes describing significant changes in each Axom release are documented
in this file.

The format of this file is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

The Axom project release numbers follow [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased] - Release date yyyy-mm-dd

### Added
- Centralized Axom's memory management functions in a seperate header and extended them to support 
  different memory spaces through Umpire. The caller may now specify the desired memory space by
  an optional argument or set a different default memory space to use. 
- Added the ability to point Axom to an UMPIRE build by specifying UMPIRE_DIR either in
  a host-config or at the command line. Axom components can link to Umpire, by specifying
  "umpire" as a dependency. A simple umpire smoke test is also added for regression testing.
- Added for_all_faces to the mint execution model.
- Added support for face connectivity in the mint UnstructuredMesh class.
- Added support for face data and face connectivity in the mint StructuredMesh classes.
- Added Python module for Quest signed distance interface.
  See the file src/axom/quest/interface/python/README.md for more information.
- Added capability in sidre IOManager to read files while running on a number
  of MPI ranks greater than the number of ranks that were used to create
  the files. 
- Users can now set the vertex welding threshold parameter in Quest's In/Out query.
  This was previously not exposed to the user. The default value is 1E-9.

### Removed
- Moved mint::Array to axom::Array with sidre storage in sidre::Array;
  also moved mint::IndexType to axom::IndexType.
- Replaced sidre::SidreLength with sidre::IndexType.
- Replaced usage of std::size_t in sidre with sidre::IndexType.
- Added AXOM_ENABLE_EXPORTS which enables CMAKE_ENABLE_EXPORTS to allow demangled
  axom function names in stack traces. This option is ON by default in debug builds.

### Deprecated

### Changed
- Replaced old quest C-style interface with a new quest inout API.
  Functions related to the inout point containment query are prefixed with "inout_".
  The new API has an option to set the verbosity of the inout initialization and query.
- Changed sidre::IndexType to be a 64bit signed integer.
- Changed slic::stack_trace to slic::internal::stack_trace which now attempts to
  output a demangled stack trace.

### Fixed
- Axom can once again be configured with `AXOM_ENABLE_EXAMPLES=ON` and `AXOM_ENABLE_TESTS=OFF`.
- Quest's vertex welding now works with small welding threshold values (e.g. 1E-20).
  Welding was previously broken when this value was smaller than 1E-8.
  This fix also resolved an issue with small grid spacing values in primal's RectangularLattice.

### Known Bugs


## [Version 0.3.0] - Release date 2018-08-02

### Added
- Added initial implementation of a RAJA-based mesh-aware execution model in Mint.
  The execution model provides a functional interface for generic traversals and
  various mesh traversals, e.g., looping over all the nodes or cells of the mesh.
- Added AXOM macros for lambda expressions and host/device decorators
- Added the ability to point Axom to a RAJA build by specifying `RAJA_DIR` at
  a host config or at the config-build.py. Axom components can now link RAJA
  by specifying "raja" as a dependency. A simple RAJA smoke test is also added
  for regression testing.
- Added `isStructured()` and `isUnstructured()` convenience methods to the mint::Mesh object.
- Added support for using MPI-3 on-node shared memory data-structures to store
  the input surface mesh in quest. Instead of each rank duplicating the mesh
  data (e.g. node coordinates and cell connectivity), ranks within the same
  compute node can now share and utilize the same mesh data buffer. This reduces
  the memory overhead associated with duplicating the mesh at each rank which can
  be a limiting factor when reading large meshes. This feature is currently only
  exposed in Quest's signed distance query and requires Axom to be compiled with
  MPI-3 support (i.e., setting `AXOM_USE_MPI3=ON`)
- Added the ability to call resize() on a Mint UnstructuredMesh. This functionality
  was not previously exposed in the UnstructuredMesh API.
- Added support for non-closed surface mesh input to the signed distance query.
- Added new interface for the signed distance query along with corresponding tests 
  and examples.
- Updated to [fmt version 5.1.0](https://github.com/fmtlib/fmt/releases/tag/5.1.0)
- Added AXOM_ENABLE_TESTS which is a CMake dependent option of ENABLE_TESTS
- Added AXOM_ENABLE_DOCS which is a CMake dependent option of ENABLE_DOCS
- Added AXOM_ENABLE_EXAMPLES which is a CMake dependent option of ENABLE_EXAMPLES
- Added jacobi_eigensolve() method for computing the eigenvalues and eigenvectors
  of real, symmetric matrices.
- Added matrix_norm() operator for computing matrix norms. The implementations
  supports the p1-norm, infinity-norm and frobenious matrix norm.
- Added eigen_sort() routine for sorting the supplied eigenvalues and corresponding
  eigenvectors in ascending order.
- Initial integration of Mint and Sidre. Mint can now operate on meshes 
  stored in Sidre and conform to the [computational mesh blueprint conventions](http://llnl-conduit.readthedocs.io/en/latest/blueprint.html).
- Added a sphere-sphere intersection test to Primal.
- Added a utility function to Quest to *weld* vertices in a triangle mesh that 
  are within a given tolerance. After welding, all triangles incident in a 
  vertex have the same index for that vertex. This function has been integrated
  into the ``mesh_tester`` utility.
- Added a bounded All-Nearest-Neighbor query to Quest. This query takes a list
  of point locations and associated regions, and for each point reports the 
  nearest point in a different region that is no farther than a max search 
  radius.
- Axom now exports [sparsehash version 2.0.3](https://github.com/sparsehash/sparsehash). Previously, it was only used internally.
- Added a BitSet class to Slam.
- Added a Tetrahedron primitive to Primal.
- Added an in_sphere operator to Primal, which is a predicate that
  is used extensively for Delaunay triangulations.

### Removed
- Axom no longer depends on the Boost library.
- Removed ENABLE_PYTHON CMake option. Python was only used by Shroud so restricted Python
  checks to when Shroud generation is enabled
- Removed Lua as a dependency of Axom.
- Removed signed distance query functions from the quest.hpp interface. The
  signed distance query is supported in its own exclusive interface.
- Removed AXOM_NULLPTR. Use nullptr instead.

### Deprecated
- 

### Changed
- Simplified the external constructors for the Mint UnstructuredMesh. Specifically, 
  the caller is no longer required to: (a) specify the dimension, which can be computed
  internaly based on other input arguments, and (b) specify explicitly the capacity for 
  the node coordinate and connectivity arrays. In the more common case, the capacity is
  equivalent to the associated size when constructing a mesh by supplied external buffers.
- Restructured source directory.  #includes will now mirror file structure.  For
  example, '#include "sidre/Group.hpp"' is now '#include "axom/sidre/core/Group.hpp"'.
- The root CMake file for Axom is now located in ``<axom>/src``'s root directory,
  rather than in ``<axom>``
- Prefixed all Axom CMake options with AXOM_ to avoid conflicts
- ENABLE_SPARSEHASE -> AXOM_ENABLE_SPARSEHASH
- ENABLE_ALL_COMPONENTS -> AXOM_ENABLE_COMPONENTS
- ENABLE_<component name> -> AXOM_ENABLE_<component name>
- MINT_USE_64BIT_INDEXTYPE -> AXOM_MINT_USE_64BIT_INDEXTYPE
- MINT_USE_SIDRE -> AXOM_MINT_USE_SIDRE
- CMake minimum is now 3.8 for non-CUDA builds and 3.9 for CUDA builds
- Axom now requires a C++11 compiler.
- Refactored Axom's Matrix/Vector operators and consolidated them in one file.
- Removed overloaded arithmetic operators from the Matrix class to avoid 
  potential negative performance impacts. Applications should use the new 
``matvecops`` methods for such operations.
- Quest STL reader now returns a status code, indicating whether reading
  the STL file was successful. Also, the reader now leverages the improved
  Mint API to reserve storage for mesh and avoid re-allocations when reading
  the STL mesh data into the Mint mesh.
- Refactored and cleaned up Primal's Sphere class.
- Refactored Mint and removed all STL usage in preparation for GPUs.

### Fixed
- Fixed minor memory leak in quest fortran example
- Bugfix for "multiply-defined" linker error in slam::Bitset and quest::PointInCellTraits

### Known Bugs
-


## [Version 0.2.9] - Release date 2018-03-08

### Added
- Updated to [conduit version 0.3.1](https://github.com/LLNL/conduit/tree/v0.3.1)
- Updated to [shroud version 0.8.8](https://github.com/LLNL/shroud/tree/v0.8.0)
- Improved platform support for LLNL's ``blue_os`` and ``bg/q`` systems.
  Axom now builds with Fortran enabled using the xlc and clang compilers.
- Improved support for Axom on Windows, including new host-configs for
  Microsoft's Visual Studio compiler and for the intel compiler on Windows.
  All Axom components can now be built on Windows, but we do not yet support
  hdf5 or Fortran on Windows.
- Added geometric Plane primitive to Primal. The Plane defines an oriented
  plane in 2D and 3D and provides support for operations such as projection of
  a point to a plane, signed distance and orientation.
- Added ability to configure Axom (in particular Sidre and Spio) without hdf5.
- Improved testing of [Scalable Checkpoint Restart (SCR)] library in Sidre.
- Added a Point-In-Cell query to Quest. The Point In Cell query finds the cell
  in a computational mesh that contains an arbitrary point in space.
  If such a cell exists, it also finds the isoparametric coordinates of the
  point with respect to the cell. The query supports higher order
  [mfem](http://mfem.org) meshes.
- Added cross-product and linspace operators to the vector utilities in 
``numerics``

### Removed
- 

### Deprecated
- 

### Changed
- The root cmake file for Axom is now located in ``<axom>``'s root directory,
  rather than in ``<axom>/src``
- ``primal`` is no longer a header-only library.
- Modified ``quest`` API to allow using a ``mint`` mesh that is already
  resident in memory.

### Fixed
- Fixed a divide-by-zero problem in ``primal::intersect()``
- Fixed the calculation of the Jacobian in ``mint::FiniteElement`` to support
  elements that are in higher-dimensional ambient space, e.g., surface elements,
  a Triangle or Quad in 3D.

### Known Bugs
-


[Unreleased]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.9&sourceBranch=refs%2Fheads%2Fdevelop&targetRepoId=1066
[Version 0.2.9]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.8&sourceBranch=refs%2Ftags%2Fv0.2.9&targetRepoId=1066

[Scalable Checkpoint Restart (SCR)]: https://computation.llnl.gov/projects/scalable-checkpoint-restart-for-mpi
