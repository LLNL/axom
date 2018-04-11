Axom {#mainpage}
================

Axom provides libraries that address common computer science needs.  It grew from the recognition that physics codes at Lawrence Livermore National Laboratories had encountered the need for logging, a data store, geometric computation, file I/O, and other facilities and had addressed the needs in various idiosyncratic methods.  Axom is designed to provide a drop-in replacement for code that wants to adopt it.

# Axom components

* @subpage axomutiltop provides shared utility functionality to all components.
* @subpage lumberjacktop provides logging aggregation and filtering capability.
* @subpage minttop provides a mesh representation with finite element operations.
* @subpage primaltop provides an API for geometric primitives and computational geometry tests.
* @subpage questtop provides an API to query point distance and position relative to meshes.
* @subpage sidretop provides a data store with hierarchical structure.
* @subpage slamtop provides an API to construct and process meshes.
* @subpage slictop provides logging levels and targets.


# Further documentation

- [User manual](../../../sphinx/web_main_docs/html/index.html)
- [Quick start guide](../../../sphinx/quickstart_guide_docs/html/index.html)
- [Development guide](../../../sphinx/dev_guide_docs/html/index.html) (intended for Axom contributors)
- Coding [style guide](../../../sphinx/coding_guide_docs/html/index.html) (intended for Axom contributors)

Each component contains a `test` subdirectory containing numerous code unit tests and an `example` subdirectory with more extensive demonstrations of the Axom library.  The Sidre, Mint, and Slam examples include several versions of the `lulesh` mini-app,  an implementation of the heat equation and a 1-D shock tube simulation.
