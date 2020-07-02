Slam {#slamtop}
=========

Axom's [Slam](@ref axom::slam), (S)et-theoretic (L)ightweight (A)PI for (M)eshes, component provides a collection of high performance, thoroughly tested, reusable components that can be combined to define distributed mesh data structures for multiphysics simulation codes.  Slam's classes and functions provide context to a mesh's entities and associated data stored as raw data.  Slam models meshes in terms of three set-theoretic concepts:

* Sets of entities (e.g. vertices, cells, particles)
* Relations among pairs of sets (e.g. incidence, adjacency and containment relations)
* Maps defining fields and attributes on the elements of a given set


<!--    (see ['components' section](@ref #components) for more detail) -->

# Design goals {#goals}

The concepts in this component are not new. In fact, since they lie at the foundation of computational meshes, they are implemented (to some degree) in every simulation code.  However, the underlying abstractions are typically latent in the design, and developers are often wary of making changes to working code.   Explicitly modeling these underlying abstractions improves the comprehensibility and maintainability of the code and presents opportunities for optimization.
<!-- (e.g. we can define some constants at compile time, when they are known). -->

Our template-based implementation helps reduce the abstraction cost (or at least moves it to compile time).

[Upcoming]
We intend to shield our users from these details through a generative model which describes the required pieces and a separate for a configurator to generate.

# Sets {#sets}

We model the basic entities of a distributed mesh, such as its vertices and cells, as *sets*, collections of entities. For performance, we implement these as [ordered sets](@ref axom::slam::OrderedSet)
where we can associate an index with each element.

Each element of such a set can be described in terms of an offset, a stride and an indirection.
That is, for a set `theSet`, the element at position `p` is:

    someSet[p] = indirection( stride*pos + offset)

<!-- ## Subsets {#subsets} -->

# Relations {#relations}

A relation is a subset of the Cartesian product of two sets.
For simulation meshes, we are often interested in a relation operator -- a function from Set A to Set B which returns the elements from the second set that are associated with the given element *x* in the first set.  We can denote this as \f$ R_A(x) = \{ y \in B | x ~_R y \} \f$ 

A taxonomy of relations
* Storage: Implicit vs. Explicit
* Mutability: Static vs. Dynamic
* Cardinality: Fixed vs. Variable


# Maps {#maps}

Maps associate a value with a set member, storing field data.


# Expected usage
* Codes will use this layer and the components (described below) to define a concrete mesh implementation.
          
* A generic outer layer that axom components can target.
  This layer will define a generic API for accessing the mesh geometry.
  All mesh implementations must satisfy this API (either natively or through an adaptor) to work with axom components.  
          
<!--  We envision a CMI-like interface for this layer.
  [It should also be possible to have an ITAPS iMesh interface as well since the design goals will be similar at this level]. -->
  
  
  
  <!--
-- example [OrderedSet](@ref axom::slam::OrderedSet)
-->
  
