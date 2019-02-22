
Slam User Documentation
=========================

The Slam (Set-theoretic Lightweight API for Meshes) component of Axom provides
high performance reusable components for building up distributed mesh data structures 
in HPC simulation codes. 

.. Goals and target audience

Slam targets the low level implementation of the underlying pieces that comprise distributed mesh data structures. 
It is aimed at developers who implement mesh data structures within HPC applications.

Its design is based on the observation that despite the apparent differences in the high level view of mesh 
data structures, many of the core features are shared at a lower level where we need to iterate over 
the mesh's index space.

Slam provides a simple, intuitive, API centered around a set-theoretic abstraction for meshes and associated data.
Specifically, it models meshes in terms of three core set-theoretic concepts: 

* **Sets** of entities (e.g. vertices, cells, particles)
* **Relations** among a pair of sets (e.g. incidence, adjacency and containment relations)
* **Maps** defining fields and attributes on the elements of a given set. 


The goal is for users to program against this interface without having to pay attentions to the specific design choices,
such as the memory layout and data containers for the underlying mesh data. The exposed API is intended to feel natural 
to end users (e.g. application developers and code physicists) who operate on the meshes that are built up 
from slam's abstractions.


.. Policy-based design

The classes in Slam use a *Policy-based* design that orthogonally decomposes the feature space in an attempt
to combat the combinatorial explosion of possible design choices without sacrificing performance. 
This makes it easier to extend support for custom features that utilize or extend the basic interface.

.. For example, some design decisions are known at compile time within a given code. 
   Providing this information allows the compiler to better optimize the generated code.  

Current limitations
-------------------

* Slam is under active development with many features planned.
* Slam's classes are highly configurable.  However, it is currently difficult to set up the templated typedefs 
  for the underlying policy classes.  This can be simplified by the introduction of Generator classes which will 
  use enumerated strings to generate the necessary typedefs. 
* Slam does not yet support GPUs.



TODO:
-----

* High-level goals and sources of requirements and use cases - 
  spectrum of mesh data structure requirements: 
  regularity of grid, variety of element types, polynomial order, 
  different ghosting schemes
* Brief summary of concepts (Set, Relation, Map) - Abstraction over low-level concepts common to simulation meshes
* Summary of current limitations and future plans


**Contents:**

.. toctree::
   :maxdepth: 2
   
   first_example
   core_concepts
   implementation_details
   examples
   more
   
`Class documentation <../../doxygen/html/annotated.html>`_
