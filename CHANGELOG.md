# AXOM Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Tight integration of Mint and Sidre. Mint can now operate on meshes stored in
  Sidre and conform to the [computation mesh blueprint conventions](http://llnl-conduit.readthedocs.io/en/latest/).
- Added a sphere-sphere intersection test to Primal.
- Added a utility function to Quest to *weld* vertices in a triangle mesh that are within
  a given tolerance. After welding, all triangles incident in a vertex have the same
  index for that vertex. This function has been integrated into the ``mesh_tester`` utility.
- Adds a bounded All-Nearest-Neighbor query to Quest.  This query takes a list
  of point locations and associated regions, and for each point reports the nearest
  point in a different region that is no farther than a max search radius.
- Axom now exports [sparsehash version 2.0.3](https://github.com/sparsehash/sparsehash).
  Previously, it was only used internally.
- Added a BitSet class to Slam.
- Added a Tetrahedron primitive to Primal.
- Added an in_sphere operator to Primal.  in_sphere is a predicate that
  is used extensively for Delaunay triangulations.

### Changed
- Refactored and cleaned up Primal's Sphere class.
- Refactored Mint and removed all STL usage in preparation for GPUs.

### Removed
- Axom no longer depends on the Boost library.

## [0.2.9] - 2017-03-08

### Added
- Updated to [conduit version 0.3.1](https://github.com/LLNL/conduit/tree/v0.3.1)
- Updated to [shroud version 0.8.8](https://github.com/LLNL/shroud/tree/v0.8.0)
- Improved platform support for LLNL's ``blue_os`` and ``bg/q`` systems.
  Axom now builds with Fortran enabled using the xlc and clang compilers.
- Improved support for Axom on Windows, including new host-configs for
  Microsoft's Visual Studio compiler and for the intel compiler on Windows.
  All Axom components can now be built on Windows, but we do not yet support
  hdf5 or Fortran on Windows.
- A geometric Plane primitive in Primal. The Plane defines an oriented
  plane in 2D and 3D and provides support for operations such as projection of
  a point to a plane, signed distance and orientation.
- Adds ability to configure Axom (in particular Sidre and Spio) without hdf5.
- Improved testing of [Scalable Checkpoint Restart (SCR)] library in Sidre.
- Adds a Point-In-Cell query to Quest. The Point In Cell query finds the cell
  in a computational mesh that contains an arbitrary point in space.
  If such a cell exists, it also finds the isoparametric coordinates of the
  point with respect to the cell. The query supports higher order
  [mfem](http://mfem.org) meshes.
- Added cross-product and linspace operators to the vector utilities in ``numerics``

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

[Unreleased]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.9&sourceBranch=refs%2Fheads%2Fdevelop&targetRepoId=1066
[0.2.9]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.8&sourceBranch=refs%2Ftags%2Fv0.2.9&targetRepoId=1066

[Scalable Checkpoint Restart (SCR)]: https://computation.llnl.gov/projects/scalable-checkpoint-restart-for-mpi
