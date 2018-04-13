
Sidre User Documentation
=========================

The Sidre (Simulation data repository) component of Axom provides tools to
centralize data management in HPC applications: data description, allocation, access,
and so forth.  The goal of Sidre is efficient coordination and sharing of data:
across physics packages in integrated applications, and between applications
that provide capabilities such as file I/O, *in situ* visualization, and analysis.

Sidre's design grew out of experience with current LLNL
applications and requirements identified for new codes to run on future
architectures.  All of these codes must carefully manage data allocation and placement
carefully to run efficiently.  Capabilities in existing codes were typically
developed specifically for each code with little regard to sharing.  In
contrast, Sidre is designed from inception to be shared by different
applications.

Introduction
-------------

Sidre provides simple application-level semantics to describe, 
allocate/deallocate, and provide access to data. Currently supported
capabilities include:

.. Do we want to make keywords bold for each bullet point?

* Separate data description and allocation operations. This allows applications
  to describe data they need and then decide how best to place the data in memory.
* Multiple different "views" into a chunk of (shared) data. A Sidre view
  includes description semantics to define data type, number of elements,
  offset, stride, etc. Thus, a chunk of data in memory can be interpreted
  conceptually in different ways.
* Externally-owned "opaque" or described data. Sidre can accept a pointer to 
  an externally-allocated data object and provide access to it by name. When
  the external data is described to Sidre, interaction with it by the user
  is similar in many ways to data that Sidre owns. When the data is not 
  described (i.e., it is "opaque"), Sidre can provide access to the data 
  via a pointer, but the consumer of the pointer must know type information
  to do anything substantial with the data.
* Attributes, or metadata associated with a Sidre view.  This metadata is
  available to the user code to facilitate program logic and is also used
  in Axom to enable selective saving of data sets to disk.
* Tree-structured data hierarchies. Many mesh-based application codes 
  organize data into hierarchies of contexts (e.g., domains, regions, blocks, 
  mesh centerings, subsets of elements containing different 
  materials, etc.). Sidre supports hierarchical, tree-based organizations 
  in a simple, flexible way.
* APIs for C++, C, and Fortran along with mechanisms to ensure inter-language
  data consistency.

So far, Sidre development has focused on designing and building flexible and
powerful concepts to build on. The Sidre API includes five main concepts:

* **Datastore.** The main interface; contains a collection of Buffers, a
  collection of default Attributes, and a tree structure of Groups.
* **Buffer.**  Describes and holds a chunk of data in memory owned by Sidre.
* **Group.** Defines a tree structure like a filesystem, where Groups are like
  folders and Views are like files.
* **View.**   Describes a conceptual layout of data in memory.  Each View
  has a collection of its explicitly-set Attributes.
* **Attribute.**  Provides an item of metadata describing a View.

These concepts will be described in more detail in later sections.

At this point, Sidre supports simple data types such as scalars, strings, and 
(multidimensional) arrays of scalars. Fundamentally, Sidre does not preclude 
the use of more complex data structures, but does not currently support them
directly. Future Sidre development will expand core functionality to support
additional features such as:

* Mechanisms to associate data with memory spaces and transfer data between
  spaces.
* Support for more complex data types.
* Complex queries and actions involving Attributes.

Support for these enhancements and others will be added based on application
needs and use cases.


**Contents:**

.. toctree::
   :maxdepth: 2

   first_example
   core_concepts
   file_io
   parallel_io_concepts
   sidre_conduit

