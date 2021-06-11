
[comment]: # (#################################################################)
[comment]: # (Copyright 2017-2021, Lawrence Livermore National Security, LLC)
[comment]: # (and Axom Project Developers. See the top-level LICENSE file)
[comment]: # (for details.)
[comment]: #
[comment]: # (# SPDX-License-Identifier: BSD-3-Clause)
[comment]: # (#################################################################)


# Axom Software Release Notes

Notes describing significant changes in each Axom release are documented
in this file.

The format of this file is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/).

The Axom project release numbers follow [Semantic Versioning](http://semver.org/spec/v2.0.0.html).


## [Unreleased] - Release date yyyy-mm-dd

### Added
- Added new CMake option to allow users to turn off Axom created tools: `AXOM_ENABLE_TOOLS`
- Inlet can now log verification errors to a user-processable list instead of using SLIC
- SPIO: Added new mapping arrays to the automatically-generated Blueprint
  index to support new schema for multi-domain parallel meshes.
- Added support for optional third-party `c2c` ("contours to codes") library for parsing 2D spline data.
  `c2c` is currently only available for Axom configurations on LLNL platforms.
- Primal's `intersect(Ray, Segment)` can now return the parametric coordinates of the intersection
  along both the ray and the segment (when the intersection exists)
- Primal's `intersect(Segment, BoundingBox)` can now return the parametric coordinates bounding the
  portion of the segment contained within the BoundingBox (when the intersection exists)
- Generalizes Quest's `InOutOctree` class to work with 2D line segment meshes. Previously,
  it only worked with 3D triangle meshes

### Changed
- `MFEMSidreDataCollection` now reuses FESpace/QSpace objects with the same basis
- Harden configuration options for BLT tools (style, code quality, etc.) against accidentally being enabled for users.  Developers will
  always give a full path (e.g. `CLANGFORMAT_EXECUTABLE`)
- Inlet: `Writer`s are passed directly to `Inlet::write` instead of being registered
- `Inlet` objects can now be constructed without a user-provided `sidre::DataStore`
- Conduit version changed to v. 0.7.2
- Renames `AXOM_DEBUG_VAR` macro to `AXOM_UNUSED_VAR` since there were many cases where the latter
  was not the appropriate name. This macro elides warnings about unused variables

### Fixed
- Fixed Primal's `intersect(Ray, Segment)` calculation for Segments that do not have unit length


## [Version 0.5.0] - Release date 2021-05-14

### Added
- Added the MFEMSidreDataCollection class for describing [MFEM] meshes and associated fields.  This
  class was adapted from MFEM's SidreDataCollection and is enabled when Axom is built with MFEM
  *and* the `AXOM_ENABLE_MFEM_SIDRE_DATACOLLECTION` CMake option is enabled.
- Added `slic::setAbortFunction` to configure a custom callback when SLIC aborts.
- Added a `batched` option to quest's `InOutOctree` containment query example application.
  This uses a kernel to test for containment on an array of points.
  The query uses OpenMP threading, when available.
- Inlet: Added support for user-defined conversions from Inlet tables to user-defined
  types, and support for arrays of user-defined types
- Added compiler define `NOMINMAX` to `axom/config.hpp.in` to avoid problems with
  the Windows `min` and `max` macros.
- Added `cpp14` variant to Spack package to allow `Inlet::LuaReader` to be used more easily.
- Inlet: Added support for string-keyed associative arrays (dictionaries)
- Inlet: Added support for defining and retrieving functions from the input file
- Inlet: Added support for YAML and JSON input files
- Inlet: Added support for mixed-key (integer and string) associative arrays
- Inlet: Added support for deeply nested containers of structs
- Inlet: Added support for `void` and strings in Lua-defined functions
- Inlet: Added `get<std::vector<T>>` for retrieving arrays without index information
- Inlet: Added a new `Writer` for generating JSON schemas which can be used by text editors
  for autocompletion
- Inlet: SphinxWriter will now document the signature of function callbacks added to a schema
- Axom::Path - New class for performing basic path operations with user-selectable
  delimiter characters
- Inlet: Added a method to `inlet::Inlet` that retrieves the set of unexpected names
  in the input file
- Inlet: Added an option to mark `Container`s as strict, which fail verification when unexpected
  entries are present
- Added support in `MFEMSidreDataCollection` for registering `QFunctions`
  (data associated with quadrature points on a mesh)
- Inlet: The internal hierarchy of an `Inlet` object can be reconstructed from a Sidre group,
  excluding callback functions
- Added new overloaded version of method
  sidre::DataStore::generateBlueprintIndex to incorporate new MPI
  features in conduit and allow for generation of a blueprint index on
  an under-decomposed parallel mesh
- Added new method sidre::View::importArrayNode to import a
  conduit::Node holding array data directly into a sidre::View
- Added support for registering material and species sets in
  `MFEMSidreDataCollection`.  These correspond to
  [`matset`](https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#material-sets)s
  and
  [`specset`](https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#species-sets)s
  in the Mesh Blueprint

### Changed
- Converted [Uberenv] to a git submodule. We previously vendored a copy of this script.
- The Sidre Datastore no longer rewires Conduit's error handlers to SLIC by default. 
  It can be  explicitly rewired using the static
  `DataStore::setConduitSLICMessageHandlers()` method.
- Inlet: Changed `SchemaCreator` to an abstract class and added missing functions
- Inlet: Added ability to access the `Reader` class from `Inlet` and Sol Lua state
  from the `LuaReader` class
- Inlet: Switched accessor interface to match that of the STL with `operator[]` and
  `T get<T>()`
- Inlet: `std::shared_ptr<T>` has been replaced with `T&` in non-owning contexts
  and `std::unique_ptr<T>` in owning contexts - specifically, within Inlet's internal
  tree structure
- Unified core and SPIO unit tests into fewer executables to limit size of build directory
- Renamed `axom::slic::UnitTestLogger` to `axom::slic:SimpleLogger` because it's used in
  more than just unit tests.
- Inlet: Input file functions can now be of arbitrary signature subject to type and arity
  restrictions
- Exported all symbols on Windows by default when compiling a dynamic library
- Updated TPL `conduit` to version 0.6.0 released Nov 2, 2020.
- Updated built-in TPL `sparsehash` to version 2.0.4 released Aug 11, 2020.
- Inlet: Exposed `primal::Vector` in Lua for use in input-file-defined functions
- The `MFEMSidreDataCollection` will now reconstruct fields and the mesh when a
  datastore is `Load`ed in
- Inlet: Cleaned up `Table` interface to eliminate ambiguity and duplicated functionality
- Inlet: Renamed `DocWriter` to `Writer` and refactored its interface
- Inlet: Renamed `Table` to `Container`
- Inlet collections of mixed or incorrect type will now fail verification, even if they're
  not marked as required
- Required Inlet collections no longer fail Inlet verification if they are empty in the input file
- Inlet: `operator bool` for `Field` and `Container` has been replaced with more precise `isUserProvided`
  and `exists`, which also returns `true` if a default value was specified.
- Updated built-in TPL `fmt` to master branch snapshot, March 26, 2021.
- Inlet: `SphinxWriter` will now print only one element schema per container instead of
  printing the same schema for each element in the container
- Updated BLT to version 0.4.0 released 9 Apr 2021
- Updated MFEM to version 4.2 released 30 Oct 2020. Axom no longer requires MFEM to be built serially
- The macro for exporting symbols is now `AXOM_EXPORT` instead of `AXOM_API`
- Updated Conduit to v0.6.0
- Updated SCR to compatibility with v3.0rc1

### Fixed
- Updated to new BLT version that does not fail when ClangFormat returns an empty
  version string.  BLT/Axom now issues a warning and disables the `style` build
  target if version is unknown or wrong.
- Inlet: Apply lambda verifiers on generic containers to individual elements
  for consistency
- Inlet: Fixed a bug relating to nested table lookups of primitive arrays and functions
- Fixed a bug relating to deeply nested callback functions in Inlet
- Inlet: Always ignore primitive array elements that do not match the requested type
- Inlet: Empty structs/collections of structs with required sub-elements no longer fail
  verification
- Quest: Fixed a bug with InOutOctree for triangles that lie on faces of octree blocks
- Updated to use newer Conduit config directory
- Add support for legacy hdf5 cmake build system

## [Version 0.4.0] - Release date 2020-09-22

### Added
- Exposed the tolerance parameter `EPS` that is used to determine intersections between
  triangles in `primal:intersect()` as an optional final parameter.
- Added BVH spatial index option to the `mesh_tester` utility for calculating
  triangle-triangle intersection.
- Added `axom::execution_space< ExecSpace >::onDevice()` to check if execution
  space is on device.
- Added Axom macro `AXOM_SUPPRESS_HD_WARN` to silence host device compiler
  warnings.
- Added option to quest's `SignedDistance` class and C API to toggle whether
  the distance query computes the sign.
- Added a `batched` option to quest's signed distance query example application.
  This computes all distance queries on an array of points using a single call to `computeDistance`.
  The query uses OpenMP threading, when available.
- Added new component, Inlet, to assist in retrieving and storing data from
  an input deck.
- Added the ability to specify an [Umpire] allocator ID to use with the
  BVH. This allows the application to use a device allocator for the BVH and 
  avoid use of Unified Memory (UM) on the GPU, which can hinder perfomrmance, 
  or use a pool allocator to mitigate the latencies associated with allocation/deallocation.
  The allocator ID is specified as an optional argument to the BVH constructor.
- Added new CMake option, `AXOM_ENABLE_ANNOTATIONS`, to enable/disable code 
  annotations in Axom. Default is OFF.
- Added Axom annotation macros. The macros can be used to annotate functions,
  using the `AXOM_PERF_MARK_FUNCTION` macro, or at a more fine grain level,
  different sections of code can be annotated by wrapping them within an
  `AXOM_PERF_MARK_SECTION` block. As a first cut, this works with NVTX tools.
  However, the hooks are in place to add support for Caliper in the future. 
- Added a simple interface to NVTX that allows an application to set the color
  and category for NVTX ranges corresponding to annotated code in Axom. The
  application can now call `axom::nvtx:set_color()` and 
  `axom::nvtx::set_category()` to set the corresponding parameters respectively. 
  This facilitates in the performance evaluation by allowing developers to easily 
  filter out calls by category or visually by setting a different color to use
  in GUI tools, such as, NVVP and NSight.
- Added a portable floating_point_limits traits class, to return min(), max(), lowest() 
  and epsilon() of a `float` or `double` type. The functionality is equivalent to that provided by 
  std::numeric_limits, but, the code is host/device decorated accordingly such that it
  can also be called on the device.
- Added initial support for ray queries using the BVH. The caller may now supply a set of rays to 
  a BVH and the BVH will return a set of candidate BVH bins that intersect each ray.
- Added initial support for bounding box queries using the BVH. The caller may
  now supply a set of bounding boxes to a BVH and the BVH will return a set of 
  candidate BVH bins that intersect each bounding box.
- Added an `axom-config.cmake` file to axom's installation to streamline incorporating axom
  into user applications. See `<axom-install>/examples/axom` for example usages.
- Added [Sol] as a built-in TPL for fast and simple `C++` and `Lua` binding.
  Sol is automatically enabled when `LUA_DIR` is found. 
  The version of Sol used in this release is `v2.20.6`, which requires `C++14`.

### Removed
- Removed the `AXOM_ENABLE_CUB` option, since `CUB` is no lonher used directly in
  Axom code. Instead, we use `RAJA::stable_sort` with RAJA-v0.12.1 and fallback
  to `std::stable_sort` with older versions of RAJA and when the code is built
  without RAJA.

### Deprecated

### Changed
- Updated Axom to support RAJA-v0.12.1 and Umpire-v4.01, but the code remains
  backwards compatible with previous versions of RAJA and Umpire.
- Transitioned Axom's code formatting tool from `Uncrustify` to [clang-format].
  Axom's clang-format rules depend on clang 10.
- Modified the command line interface for `mesh_tester` utility. Interface
  now uses a *-m, --method* option to select the spatial index, and *-p, policy*
  option now accepts a string or integer value.
- Renamed the `AXOM_USE_MPI3`option to `AXOM_ENABLE_MPI3` for consistency.
- Modified the API for the BVH to accomodate different query types. The queries are now
  more explicitly called `BVH::findPoints()` and `BVH::findRays()`.
- Modified the API of Axom's memory management routines to not leak usage of Umpire. Instead of 
  passing an `umpire::Allocator` object to specify an allocator, we now use the corresponding
  integer ID associated with the allocator.
- All names in the C API now preserve the case of the C++ function.
  ex. `SIDRE_datastore_new` is now `SIDRE_DataStore_new`.
- Fortran API in slic module. `axom::slic::message` Level enums are changed
  from  *enum-name_enumerator* to *namespace_enumerator*.
  ex. `level_error` is now `message_error`.
- Fortran derived-type constructors are now generic functions named afer the derived type.
  `datastore_new` is now `SidreDataStore`
  `iomanager_new` is now `IOManager`

### Fixed
- Fixed a bug in `primal::intersect(Segment, BoundingBox)` and added regression tests.
- Spin's octrees can now be used with 64-bit indexes. This allows octrees 
  with up to 64 levels of resolution when using a 64-bit index type.
- Resolved issue with `AXOM_USE_64BIT_INDEXTYPE` configurations. Axom can once again
  be configured with 64-bit index types.
- Fixed a triangle-triangle intersection case in primal that produced inconsistent results
  depending on the order of the triangle's vertices.
- Fixed issue in the parallel construction of the BVH on GPUs, due to incoherent
  L1 cache that could result in some data corruption in the BVH. The code now
  calls ``__threadfence_system()`` after the parent is computed and stored back
  to global memory to ensure that the *write*  is visible to all threads. 
- Fixed issue in Mint that would cause the clang@9.0.0 compiler to segfault. The
  `mint_cell_types.cpp` test was causing a segfault in the compiler. The main
  issue triggering this compiler bug was the use of `constexpr` when defining the
  static `cell_info` array of structs. The code has been modified to use `const`
  instead.
- Fixed issue in Quest's Signed Distance query that would prevent consecutive
  calls to Quest when MPI-3 shared memory is enabled due to not properly 
  nullifying internal pointers when finalize is called.
- Fixed issue where the BVH would dispatch to the CPU sort() routine when the
  specified execution policy was CUDA_EXEC async. Now, when the execution policy
  is CUDA_EXEC the code would correctly dispatch to the GPU sort, using CUB
  (when CUB is enabled), regardless of whether it's synchronous or asynchronous.
- Fixed issue with missing the bvh_traverse.hpp from the install prefix, which was preventing
  applications from using the BVH when pointing to an Axom install prefix.
- Fixed usage of cuda kernel policies in Mint. Raja v0.11.0 changed the way max threads
  launch bounds is calculated. Consequently, a large number of threads was being launched
  leading to max registry count violation when linking. We are now using fixed kernel size
  of 256 threads (16x16 in 2D and 8x8x4 in 3D).
- Third-party libraries can now build on the Windows platform through uberenv using vcpkg
  ("zero-to-axom support on Windows")

### Known Bugs
- Encountered a compiler bug on IBM LC platforms when using the IBM XL C/C++
  compiler. The issue is manifested in the `generate_aabbs_and_centroids` method
  in the `spin_bvh.cpp` unit test. It seems that the compiler does not handle
  the lambda capture of the arrays correctly which leads to a segfault. A
  workaround for the IBM XL compiler is provided.
- There is a known bug in MVAPICH that prevents consecutive creation/deletion
  of MPI windows. This was encountered on LC platforms when enabling shared
  memory in the Signed Distance Query. See the corresponding 
  [Github Issue](https://github.com/LLNL/axom/issues/257) for details.

## [Version 0.3.3] - Release date 2020-01-31

### Added
- Define different execution spaces. This refines and consolidates
  the execution policy concepts from mint and spin, which are now defined in
  Axom core, such that they can be used by other components.
- Added a generic axom::for_all(), which can be used to write simple parallel
  loops.
- Added [CLI11](https://github.com/CLIUtils/CLI11) command line parser as a built-in third party library.

### Removed
  
### Deprecated

### Changed
- Updated Conduit to v0.5.1
- Updated RAJA to v0.11.0
- Updated Umpire to v2.1.0, which, natively handles zero byte re-allocations consistently. Workaround
  for older releases of Umpire is in place for backwards compatibility. 
- Updated BLT to develop (f0ab9a4) as of Jan 15, 2020
- Set CUDA_SEPARABLE_COMPILATION globally instead of just in a few components.
- Reduced verbosity of quest's InOutOctree in cases where query point lies on surface.
- Changed semantics of ``axom::reallocate(_, 0)`` to always return a valid pointer.
  This matches the semantics of Umpire's ``reallocate(_, 0)``.
  Note: Umpire's PR-292 fixed a bug in its handling of this case and Axom
  includes a workaround to get the new behavior until we upgrade to Umpire v2.0+.

### Fixed
- Fixed a bug in ``convert_sidre_protocol`` example. Data truncation functionality now
  works properly when multiple Views point to the same data.

### Known Bugs


## [Version 0.3.2] - Release date 2019-09-22

### Added
- Added support in Mint for reading and writing an unstructured mesh in the [SU2 Mesh file format].
  This includes support for both single and mixed cell type topology unstructured mesh types.
- Added a new option to enable/disable use of CUB, AXOM_USE_CUB, which is disabled by default. This
  allows to disable CUB to circumvent issues encountered with the device linker.
- Added a BezierCurve primitive type to primal. A new ``intersect`` operator was also added to
  compute the intersection points between a pair of Bezier curves of arbitrary order.

### Removed

### Deprecated

### Changed
- Updated Raja TPL to v0.9.0
- Updated Umpire TPL to v1.0.0
- AXOM_USE_OPENMP is now being set at configure time accordingly instead of
  auto-detected based on whether `_OPENMP` is passed by the compiler. This
  fixes issues where a host code would compile Axom w/out OpenMP, but, use
  Axom in parts of the code where OpenMP is enabled.

### Fixed
- Fixed usage of Umpire's MemoryResourceType enum in Axom. Axom was assuming that
  there was a one-to-one correspondance of the entries in the MemoryResourceType enum
  and the IDs of the predefined allocators. However, this assumption generally does
  not hold. This version corrects this by explicitly querying the ID of the predefined
  allocator for a particular resource and using that subsequently in the code.

### Known Bugs

## [Version 0.3.1] - Release date 2019-07-22

### Added
- Added a new implementation of the Bounding Volume Hierarchy(BVH) spatial
  acceleration data-structure to spin. The new BVH implementation relies on RAJA
  and allows for constructing and using the BVH sequentially or in parallel on
  the CPU or in parallel on the GPU. The new implementation is available only
  when Axom is compiled with RAJA and Umpire support.
- Added Umpire and RAJA to the spack build.
- Centralized Axom's memory management functions in a separate header and extended them
  to use Umpire when enabled.
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
- Unify all Axom component libraries into one library named axom.
- The routine that checks if a surface mesh is watertight now marks boundary
  cells. A cell-centered field, named "boundary", is used to mark  boundary cells
  with a value  of  "1" and "0" otherwise. This facilitates in visually inspecting the
  surface mesh and identify the problematic regions for the In/Out and SignedDistance
  queries.

### Removed
- Moved mint::Array to axom::Array with sidre storage in sidre::Array;
  also moved mint::IndexType to axom::IndexType.
- Replaced sidre::SidreLength with sidre::IndexType.
- Replaced usage of std::size_t in sidre with sidre::IndexType.
- Added AXOM_ENABLE_EXPORTS which enables CMAKE_ENABLE_EXPORTS to allow demangled
  axom function names in stack traces. This option is ON by default in debug builds.

### Deprecated

### Changed
- Updated conduit TPL to v0.4.0
- Updated mfem TPL to v4.0
- Slam's Set and Relation classes are now parameterized by a PositionType and ElementType.
  Its Map classes are now parametrized by a SetType.
- Updated the fmt tpl.
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

[Unreleased]:    https://github.com/LLNL/axom/compare/v0.5.0...develop
[Version 0.5.0]: https://github.com/LLNL/axom/compare/v0.4.0...v0.5.0
[Version 0.4.0]: https://github.com/LLNL/axom/compare/v0.3.3...v0.4.0
[Version 0.3.3]: https://github.com/LLNL/axom/compare/v0.3.2...v0.3.3
[Version 0.3.2]: https://github.com/LLNL/axom/compare/v0.3.1...v0.3.2
[Version 0.3.1]: https://github.com/LLNL/axom/compare/v0.3.0...v0.3.1
[Version 0.3.0]: https://github.com/LLNL/axom/compare/v0.2.9...v0.3.0
[Version 0.2.9]: https://github.com/LLNL/axom/compare/v0.2.8...v0.2.9

[clang-format]: https://releases.llvm.org/10.0.0/tools/clang/docs/ClangFormatStyleOptions.html
[MFEM]: https://mfem.org
[Scalable Checkpoint Restart (SCR)]: https://computation.llnl.gov/projects/scalable-checkpoint-restart-for-mpi
[SU2 Mesh file format]: https://su2code.github.io/docs/Mesh-File/
[Umpire]: https://github.com/LLNL/Umpire
[Sol]: https://github.com/ThePhD/sol2
[Uberenv]: https://github.com/LLNL/uberenv
