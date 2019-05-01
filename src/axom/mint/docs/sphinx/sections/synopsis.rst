.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level COPYRIGHT file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

Synopsis
---------

Mint is a component of the `Axom Toolkit`_, being developed at
`Lawrence Livermore National Laboratory (LLNL) <http://www.llnl.gov>`_,
consisting of a *comprehensive*, *flexible* and *extensible*
:ref:`sections/mesh_data_model`, implemented in C++.
Mint aims to facilitate the development of large-scale, massively parallel,
multi-physics applications. Towards that goal, Mint serves as a
foundational building block that underpins the development of computational
tools and numerical discretization methods, enabling implementations that are
born *parallel* and *portable* to new and emerging architectures.

.. topic:: Salient Features

  Some of the salient features of Mint include, but are not limited to, the
  following:

  * Support for 1D/2D/3D mesh geometry.

  * Efficient data structures to represent :ref:`ParticleMesh`,
    :ref:`StructuredMesh` and :ref:`UnstructuredMesh` types, including
    unstructured meshes with :ref:`MixedCellTopology`.

  * Native support for a variety of commonly employed :ref:`CellTypes`.

  * A flexible :ref:`MeshStorageManagement` system, which can optionally
    inter-operate with `Sidre`_ as the underlying, in-memory, hierarchical
    datastore, facilitating the integration across packages.

  * Basic support for :ref:`sections/fem`, consisting of
    commonly employed :ref:`ShapeFunctions` and :ref:`Quadratures`.

  * A Mesh-Aware :ref:`sections/execution_model`, based on the `RAJA`_ programming
    model abstraction layer, enabling the implementation of computational kernels
    that are born parallel and portable to emerging architectures.

.. topic:: Requirements



  Mint is designed to be *light-weight* and *self-contained*.
  The only requirement for using Mint is a C++11 compliant compiler.
  However, to realize the full spectrum of capabilities provided,
  the following third-party libraries are supported:

  * `RAJA`_, used for the parallel execution and portability layer.

  * `Conduit`_ , for using `Sidre`_ as the :ref:`MeshStorageManagement` system.

  * `Umpire`_, for memory management on next-generation architectures.

  For further information on how to build the `Axom Toolkit`_ using these
  third-party libraries, consult the `Axom Quick Start Guide`_

.. topic:: About this Guide


  This guide discusses the basic concepts and overarching design architecture of
  Mint. The :ref:`sections/quick_introduction` section offers an excellent resource
  for getting started with Mint quickly, illustrating basic concepts
  and key capabilities in the context of a working code example.
  Additional examples and code snippets are provided in the
  :ref:`sections/tutorial` section that are designed to demonstrate specific Mint
  concepts and capabilities in a structured and simple format. For more detailed
  documentation of the interfaces of the various classes and functions in Mint,
  developers are advised to consult the `Mint Doxygen API Documentation`_. Last,
  complete examples and code walk-throughs of mini-apps using Mint are provided in
  the :ref:`sections/examples` section.

  For any questions, please consult the :ref:`sections/faq` section or feel free
  to post your question to the Axom Developers mailing list at axom-dev@llnl.gov.

