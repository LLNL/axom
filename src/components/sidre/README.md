Sidre: Simulation Data Repository {#sidretop}
=========

The ASC toolkit provides Sidre, a data repository library for flexible, low-overhead organization of buffers in a hierarchical tree.  Four concepts are embodied in Sidre classes:

* Buffers are linear arrays of data.
* Views provide access to data in buffers, specifying type, offset, and stride.
* Groups contain views.  Groups also contain other groups, in a tree structure.
* Stores contain buffers and a root group.

# Work flow {#workflow}

A code will instantiate a Sidre data store and create a hierarchy of groups, then allocate buffers and create views (within the group hierarchy) into the buffers.  Each view specifies the type and shape of the data, as well as offset and stride, and provides the mechanism for access to the buffers' contents.  A view can also provide access into a memory buffer allocated by an outside code.  Alternately, the SPIO component of the ASC toolkit can produce a Sidre store from an external data source.

After instantiantion, a code can use a string, formatted for proper path notation, to obtain a pointer to a group or view.  A view provides access to its underlying buffer, 

# Design goals {#goals}

Sidre was designed to provide a data store library that is conceptually close to existing functionality in current physics codes, and easily mapped to storage options such as HDF5.  Sidre also aims to be simple.  In contrast to HDF5, all buffers are contained by a store directly, not by groups, and groups are arranged in a tree structure, not a directed acyclic graph.  This simplicity promotes robust code and ease of use, while remaining powerful enough to fulfill the needs of physics codes.

