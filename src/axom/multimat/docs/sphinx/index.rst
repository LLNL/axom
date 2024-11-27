.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Multimat User Guide
====================

Axom's MultiMat component is a data management library for multimaterial field data
within multiphysics simulation codes. Simulation codes use materials to add extra
parts and details into the mesh without requiring those features to be modeled
conformally by adding zones that are devoted to the part geometry. This method is
often called "shaping" and is covered in Axom's Klee and Quest components. 

Multimat can be used to designate a set of materials and how much of each material
is present in each cell of a mesh. Cells can contain multiple materials, each with
an associated volume fraction that indicates how much of a given cell is occupied
by the material. 

In addition to representing materials on a mesh, MultiMat can use its material knowledge
to define fields on the mesh, and over material subsets of the mesh where multiple
values are needed when a cell contains multiple materials. MultiMat supports flexible
data mappings, layouts, and dense vs sparse field storage, allowing field data to
occupy less memory than would be possible using dense field storage.



API Documentation
-----------------

Doxygen generated API documentation can be found here: `API documentation <../../../../doxygen/html/coretop.html>`_


.. toctree::
   :caption: Contents
   :maxdepth: 1

   multimat_materials
   multimat_field_concepts
   multimat_using_fields
