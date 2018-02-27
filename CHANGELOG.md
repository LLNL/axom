# AXOM Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- A Plane geometric primitive is added in Primal. The Plane defines an oriented
  plane in 2D and 3D and provides support for operations such as, projection of
  a point to a plane, signed distance and orientation. 
- Adds ability to configure Axom (in particular Sidre and Spio) without hdf5. 
- Adds a Point-In-Cell query to Quest. The Point In Cell query finds the cell
  in a computational mesh that contains an arbitrary point in space.
  If such a cell exists, it also finds the isoparametric coordinates of the 
  point with respect to the cell. The query supports higher order
  mfem meshes (mfem.org).


[Unreleased]: https://lc.llnl.gov/bitbucket/projects/ATK/repos/axom/compare/commits?targetBranch=refs%2Ftags%2Fv0.2.8&sourceBranch=refs%2Fheads%2Fdevelop&targetRepoId=1066
