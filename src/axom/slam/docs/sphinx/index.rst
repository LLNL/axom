
Slam User Documentation
=========================

The Slam (Set-theoretic Lightweight API for Meshes) component of the CS Toolkit provides
a collection of high performance, thoroughly tested, reusable components 
that can be combined to define distributed mesh data structures 
for multiphysics simulation codes. 

The Slam component is a collection of low-level primitives 
that provide context to the mesh's entities and associated data stored as raw data. 
Specifically, it models meshes in terms of three set-theoretic concepts: 

* [Sets](#sets) of entities (e.g. vertices, cells, particles)
* [Relations](#relations) among pairs of sets (e.g. incidence, adjacency and containment relations)
* [Maps](#maps) defining fields and attributes on the elements of a given set. 



Introduction
------------

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
   
`Class documentation <../../doxygen/html/annotated.html>`_
