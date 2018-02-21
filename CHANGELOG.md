# AXOM Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.9] - 2017-02-21

### Added
- A Plane geometric primitive is added in Primal. The Plane defines an oriented
  plane in 2D and 3D and provides support for operations such as, projection of
  a point to a plane, signed distance and orientation. 
- Adds ability to configure Axom (in particular Sidre and Spio) without hdf5. 
- Improves support for building Axom on Windows, including new host-configs for
  Microsoft's Visual Studio compiler and for the intel compiler on Windows.
  All Axom components can now be built on Windows, but we do not yet support
  hdf5 or fortran on Windows.
- Improved testing of [Scalable Checkpoint Restart (SCR)] library in Sidre.
- Adds a Point-In-Cell query to Quest. The Point In Cell query finds the cell
  in a computational mesh that contains an arbitrary point in space.
  If such a cell exists, it also finds the isoparametric coordinates of the 
  point with respect to the cell. The query supports higher order
  mfem meshes (mfem.org).
  
### Changed
- The root cmake file for Axom is now located in ``<axom>``'s root directory, 
  rather than in ``<axom>/src``
- ``primal`` is no longer a header-only library.
- Modified ``quest`` API to allow using a ``mint`` mesh that is already 
  resident in memory.

### Fixed
- Fixed a divide-by-zero problem in ``primal::intersect()``

[Unreleased]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.8&sourceBranch=refs%2Fheads%2Fdevelop&targetRepoId=1066
[0.2.9]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.8&sourceBranch=refs%2Fheads%2Fdevelop&targetRepoId=1066
[Scalable Checkpoint Restart (SCR)]: https://computation.llnl.gov/projects/scalable-checkpoint-restart-for-mpi
