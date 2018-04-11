Sidre: Simulation Data Repository {#sidretop}
=========

Axom provides Sidre, a data repository library for flexible, low-overhead organization of buffers in a hierarchical tree.  Four concepts are embodied in Sidre classes:

* [Buffers](@ref axom::sidre::Buffer) are linear arrays of data.
* [Views](@ref axom::sidre::View) provide access to data in Buffers, specifying type, offset, and stride.
* [Groups](@ref axom::sidre::Group) contain Views.  Groups also contain other Groups, in a tree structure.
* [Data stores](@ref axom::sidre::DataStore) contain Buffers and a root group.

Datastores also maintain a list of [Attributes](@ref axom::sidre::Attribute) to hold metadata to apply to Views.

# Work flow {#workflow}

A code will instantiate a Sidre data store and create a hierarchy of groups, then allocate buffers and create views (within the group hierarchy) into the buffers.  Each view specifies the type and shape of the data, as well as offset and stride, and provides the mechanism for access to the buffers' contents.  A view can also provide access into a memory buffer allocated by an outside code.  Alternately, a Sidre store can be read in from an external data source.  After instantiantion, a code can obtain a pointer to a group or view by path name.

The [Sidre guide](../../../sphinx/sidre_docs/html/index.html)
introduces these concepts in more detail.

