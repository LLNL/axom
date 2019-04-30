
Primal User Documentation
=========================

Primal is a component of Axom that provides efficient and general purpose
algorithms and data structures for computational geometry.  Primal provides:

- Classes to represent geometric primitives such as Point and Ray
- Functions operating on Primal's classes to implement geometric operators,
  including distance and intersection
- Two spatial query acceleration data structures used to accelerate operations
  such as collision detection and containment.

This tutorial contains a collection of brief examples demonstrating
Primal primitives and operators.  The examples instantiate geometric primitives
as needed, demonstrate geometric operators, and build the UniformGrid 
and BVHTree spatial data structures.  These examples show representative 
overloads of each of the Primal operators (see the
`API documentation <../../../../doxygen/axom_doxygen/html/primaltop.html>`_ 
for more details).

.. toctree::
   :maxdepth: 2

   primitive
   operator
   spatialindex
