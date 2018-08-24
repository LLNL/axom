
Primal User Documentation
=========================

Primal is a component of Axom that provides efficient implementations of 
fundamental geometric ideas.  Primal provides:

- Classes to represent geometric primitives such as Point and Ray
- Functions on Primal's classes implementing geometric operators, including
  distance and intersection
- Classes implementing two spatial indexes, data structures used to accelerate operations
  such as collision detection and containment.

Here is a collection of introductory examples showing Primal primitives and
operators.  The examples instantiate geometric primitives as needed, demonstrate
geometric operators, and build the UniformGrid and BVHTree spatial index
objects around them.  These examples show representative overloads of each of
the Primal operations (see the 
`API documentation <../../../doxygen/axom_doxygen/html/primaltop.html>`_ 
for more details).

.. toctree::
   :maxdepth: 2

   primitive
   operator
   spatialindex
