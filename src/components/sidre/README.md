Sidre: Simulation Data Repository {#sidretop}
=========

Axom provides Sidre, a data repository library for flexible, low-overhead organization of buffers in a hierarchical tree.  Five concepts are embodied in Sidre classes:

* [Buffers](@ref axom::sidre::Buffer) are linear arrays of data.
* [Views](@ref axom::sidre::View) provide access to data in Buffers, specifying type, offset, and stride.
* [Groups](@ref axom::sidre::Group) contain Views.  Groups also contain other Groups, in a tree structure.
* [Data stores](@ref axom::sidre::DataStore) contain Buffers and a root group.
* [Attributes](@ref axom::sidre::Attribute) hold metadata that apply to Views.

# Work flow {#workflow}

A code will typically
- Instantiate a Sidre data store
  - create a hierarchy of groups
  - allocate buffers and create views into the buffers
  - create a list of attributes and set attribute values on views
- Alternately, read in a data store from an external data source
- Retrieve groups and views from the data store by name
- Using the pointer from a view, read from and write into data arrays
- Save the hierarchical data store to a file, such as an HDF5 archive
  - Save the subset of the data store: only the views with a specific attribute set
- Dispose of the data store
  - External pointers must be dealt with manually

The [Sidre guide](../../../sphinx/sidre_docs/html/index.html)
introduces these concepts in more detail.

