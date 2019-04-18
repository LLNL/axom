.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)


.. _sections/execution_model:

Execution Model
================

Mint provides a mesh-aware :ref:`sections/execution_model`, that
is based on the `RAJA`_ programming model abstraction layer, enabling the
implementation of computational kernels that are born *parallel* and *portable*
to emerging architectures.

.. note::
   To utilize NVIDIA GPUs, using the `RAJA`_ CUDA backend, Axom needs to be
   compiled with CUDA support and linked to a CUDA-enabled `RAJA`_ library.
   Consult the `Axom Quick Start Guide`_ for more information.

The execution model consists of a set of templated functions that accept a
mesh instance and a C++11 Lambda Expression as an argument. The Lamda
expression encapsulates the body of the computational kernel that operates
on the supplied mesh. The general form of the constituent templated
functions of the :ref:`sections/execution_model` is shown in
:numref:`figs/execModel`.

.. _figs/execModel:
.. figure:: ../figures/execmodel.png
   :align: center
   :scale: 50%
   :alt: Execution Model

   General form of the constituent templated functions of the
   :ref:`sections/execution_model`

As shown in :numref:`figs/execModel`, the key elements of the functions
that comprise the :ref:`sections/execution_model` are:

* **Iteration Space:** the iteration space is defined by the suffix of the
  function and defines the mesh entity that is being traversed, e.g., the
  :ref:`Nodes`, :ref:`Cells` or :ref:`Faces` of the mesh.

* **Execution Policy:** the execution policy is defined as a template argument
  and it specifies *where* and *how* the kernel is executed. Mint defines a
  set of high-level :ref:`executionPolicies` that map to corresponding
  `RAJA`_ execution policies internally.

* **Execution Signature:** the execution signature is an optional template
  argument and it specifies any additional arguments that the supplied
  kernel, which is encapsulated in a C++ Lambda expression takes.

* **Kernel:** the kernel is defined at the application layer, by a C++11
  Lambda expression, and is given as an argument to the function of the
  :ref:`sections/execution_model`

.. note::
    By default, if a second template argument is not specified to any of the
    constituent functions of the :ref:`sections/execution_model`, the kernel,
    usually specified by a `Lambda Expression`_, is assumed to take a single
    argument that corresponds to the iteration index.


.. _executionPolicies:

Execution Policies
------------------

.. _traversals:

Traversals
------------

.. _genericTraversals:

Generic Traversals
^^^^^^^^^^^^^^^^^^

.. _NodeTraversals:

Node Traversals
^^^^^^^^^^^^^^^

.. code-block:: cpp

   mint::Mesh* domain = get_mesh();
   double* phi = domain->getFieldPtr< double >( "phi", mint::NODE_CENTERED );

   ...

   mint::for_all_nodes< exec >( domain, AXOM_LAMBDA(IndexType nodeIdx) {

     ...

   } );


.. _CellTraversals:

Cell Traversals
^^^^^^^^^^^^^^^

.. code-block:: cpp

   mint::Mesh* domain = get_mesh();
   double* phi = domain->getFieldPtr< double >( "phi", mint::CELL_CENTERED )

   ...

   mint::for_all_cells< exec >( domain, AXOM_LAMBDA(IndexType cellIdx) {

     ...

   } );


.. _FaceTraversals:

Face Traversals
^^^^^^^^^^^^^^^

.. code-block:: cpp

   mint::Mesh* domain = get_mesh();
   double* phi = domain->getFieldPtr< double >( "phi", mint::FACE_CENTERED )

   ...

   mint::for_all_faces< exec >( domain, AXOM_LAMBDA(IndexType faceIdx) {

     ...

   } );


.. warning::
    This section is under development.

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
