Sidre {#sidretop}
=========

The Sidre component in Axom provides a data repository library for flexible, low-overhead, hierarchical organization of simulation data. Five concepts are embodied in Sidre classes:

* [Buffer](@ref axom::sidre::Buffer) is a container for data in memory.
* [View](@ref axom::sidre::View) provides a virtual description of data referenced through a pointer (type, offset, stride) and access to that data.
* [Group](@ref axom::sidre::Group) is a node in a hierarchical tree structure for data. A group may contain any number of (child) Groups or Views.
* [DataStore](@ref axom::sidre::DataStore) is the central access point for a Sidre data hierarchy. It contains Buffers and a root Group.
* [Attribute](@ref axom::sidre::Attribute) holds metadata that apply to Views for selective processing of data.

Sidre also provides parallel I/O facilities via the [IOManager](@ref axom::sidre::IOManager) class.

# Work flow {#workflow}

Typically, an application code will
- Create a DataStore object
  - Create a hierarchy of groups and views in the data store; the group and view hierarchy, once created, can be modified as needed by adding, destroying, moving, or copying groups and views within it
  - Describe data associated with views
  - Allocate data in buffers or in views directly (there are many ways to do this)
  - Create a set of attributes and associate attribute (and values) with views
- Alternately, read in a data store hierarchy and data from an external source, such as a set of HDF5 files
- Access groups and views in a hierarchy by name (path-like syntax is supported to access any group or view in a subtree of any group in the hierarchy)
- Access a data pointer from a view to read or write data in memory
  - Pointers to external data allocations can also be managed
- Save data in a group hierarchy files, such as an HDF5 archive
  - Save a subset of data in a group hierarchy; e.g., only the views with a specific attribute set
- Delete the data store when done
  - This will also clean up data allocated via Sidre mechanisms; external data must be deallocated by user code
