Axom {#mainpage}
================

Axom provides libraries that address common computer science needs.  It grew from the recognition that physics codes at Lawrence Livermore National Laboratories had encountered the need for logging, a data store, geometric computation, file I/O, and other facilities and had addressed the needs in various idiosyncratic methods.  Axom is designed to provide a drop-in replacement for code that wants to adopt it.

# Axom components

* @subpage coretop provides shared utility functionality to all components.
* @subpage lumberjacktop provides logging aggregation and filtering capability.
* @subpage minttop provides a mesh representation with finite element operations.
* @subpage primaltop provides an API for geometric primitives and computational geometry tests.
* @subpage questtop provides an API to query point distance and position relative to meshes.
* @subpage sidretop provides a data store with hierarchical structure.
* @subpage slamtop provides an API to construct and process meshes.
* @subpage slictop provides logging levels and targets.
* @subpage spintop provides spatial acceleration data structures, also known as spatial indexes.

Dependencies between components are as follows:
- Core has no dependencies, and the other modules depend on Core
- Slic optionally depends on Lumberjack
- Slam, Primal, Mint, Quest, Spin, and Sidre depend on Slic
- Mint optionally depends on Sidre
- Spin depends on Primal and Slam
- Quest depends on Slam, Primal, Spin, and Mint

The figure below summarizes the dependencies between the modules.  Solid links
indicate hard dependencies; dashed links indicate optional dependencies.

\dotfile dependencies.dot "Module dependencies" width=430px


# Further documentation

- [User manual](../../../sphinx/web_main_docs/html/index.html)
- [Quick start guide](../../../sphinx/quickstart_guide_docs/html/index.html)
- [Development guide](../../../sphinx/dev_guide_docs/html/index.html) (intended for Axom contributors)
- Coding [style guide](../../../sphinx/coding_guide_docs/html/index.html) (intended for Axom contributors)
- Please see our [license](../../../LICENSE) and accompanying [notice](../../../NOTICE).

Each component contains a `test` subdirectory containing numerous code unit tests and an `example` subdirectory with more extensive demonstrations of the Axom library.  The Sidre, Mint, and Slam examples include several versions of the `lulesh` mini-app,  an implementation of the heat equation and a 1-D shock tube simulation.
