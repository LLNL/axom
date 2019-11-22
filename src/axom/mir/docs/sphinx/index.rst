.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

=======================
Mir User Documentation
=======================

Axom's Material Interface Reconstruction (MIR) component provides algorithms for
reconstructing the interface surfaces between different materials in multimaterial
meshes. The goal of the component is to provide simulation code developers with a 
set of algorithms for performing material interface reconstruction with efficient 
implementations and an easy-to-use API.

This component currently provides support for two MIR algorithms:

.. toctree::
   :maxdepth: 2

   equi_z
   iterative_equi_z

.. raw:: html

    <h3>Current limitations</h3>

.. note:: This MIR component is under active development with many features planned.

* Improved efficiency for the Equi-Z algorithm implementation is under development.
* Support for GPUs in MIR is planned.
* Support for material interface reconstruction in high-order meshes is planned.


.. raw:: html

    <h3>Example Code</h3>

In order to access the MIR functionality in this component, import the MIR header file:

.. literalinclude:: ../../examples/mir_minimal_tutorial.cpp
   :start-after: _mir_header_start
   :end-before: _mir_header_end
   :language: C++

The following is all the code needed to generate a mesh with the ``MIRMesh`` class from one of 
the default test cases, perform material interface reconstruction using the Equi-Z algorithm 
with the ``computeReconstructedInterface()`` function, and then write out the resulting mesh 
to a vtk file.

.. literalinclude:: ../../examples/mir_minimal_tutorial.cpp
   :start-after: _mir_main_loop_start
   :end-before: _mir_main_loop_end
   :language: C++