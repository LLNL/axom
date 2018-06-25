.. ##
.. ## Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

.. _sections/mesh_data_model:

.. _MeshDataModel:

Mesh Data Model
===============

This section presents the key constituents of Mint's :ref:`MeshDataModel`.
First, :ref:`PreliminaryConcepts` and a high-level description of the underlying
:ref:`MeshRepresentation` in Mint are presented, followed by,
a classification and overview of the different :ref:`MeshTypes`.
These concepts are then combined in presenting the :ref:`Architecture` of
the Mint :ref:`MeshDataModel` and underlying :ref:`MeshStorageManagement`
substrate.

.. _PreliminaryConcepts:

Preliminary Concepts
---------------------

A mesh [#f1]_, denoted by :math:`\mathcal{M}(\Omega)`, provides a discrete
represenation of a geometric domain of interest, :math:`\Omega`, on which, the
underlying *mathematical model* is evaluated. The mathematical model
is typically defined by a system of governing *Partial Differential Equations
(PDEs)* and associated boundary and initial conditions. The solution
to the governing PDE predicts a physical process that occurs and evolves on
:math:`\Omega` over time. For example, consider the flow around an aircraft,
turbulence modeling, blast wave propagation over complex terrains, or,
heat transfer in contacting objects, to name a few.
Evolving the mathematical model to predict such a physical process is typically
done numerically, which requires discretizing the governing PDE by a numerical
scheme, such as, a Finite Difference (FD), Finite Volume (FV), or, the
Finite Element Method (FEM), chief among them.

.. _figs/meshedDomain:
.. figure:: ../figures/meshed_domain.png
   :align: center
   :scale: 100%
   :alt: Sample Mesh domain

   Mesh discretization of a geometric domain: (a) Sample geometric domain,
   :math:`\Omega`. (b) Corresponding mesh of the domain,
   :math:`\mathcal{M}(\Omega)`. The *nodes* and *cells* of the mesh, depicted in
   red, correspond to the discrete locations where the unknown variables of the
   governing PDE are stored and evaluated.

Discretization of the governing PDE requires the domain to be approximated
with a mesh. For example, :numref:`figs/meshedDomain` (a) depicts a geometric
domain, :math:`\Omega`. The corresponding mesh, :math:`\mathcal{M}(\Omega)`,
is illustrated in :numref:`figs/meshedDomain` (b). The mesh approximates
the geometric domain, :math:`\Omega`, by a finite number of simple geometric
entities, such as, *nodes* and *cells*, depicted in red in
:numref:`figs/meshedDomain` (b). These geometric entities comprising the mesh
define the discrete locations, in space and time, at which the unknown variables,
i.e., the *degrees of freedom* of the governing PDE, are evaluated, by the
numerical scheme being employed.

There are a variety of different :ref:`MeshTypes` one can choose from.
The type of mesh employed depends on the choice of the underlying
numerical discretization scheme. For example, a finite difference scheme
typically requires a :ref:`StructuredMesh`. However, the finite volume and
finite element methods may be implemented for both :ref:`StructuredMesh` and
:ref:`UnstructuredMesh` types. In contrast, *meshless* or *mesh-free* methods,
such as, *Smoothed Particle Hydrodynamics (SPH)*, discretize the governing PDE
over a set of *particles* or *nodes*, using a :ref:`ParticleMesh`
representation.

.. #############################################################################
..  MESH Representation
.. #############################################################################

.. _MeshRepresentation:

Mesh Representation
--------------------

Irrespective of the mesh type, a mesh essentially provides a data structure
that enables efficient *storage*, *management* and *access* of:

#. Mesh :ref:`Geometry` information (i.e., nodal coordinates),

#. Mesh :ref:`Topology` information (i.e., cell-to-node connectivity, etc.), and

#. :ref:`FieldData` stored on a mesh

The underlying, concrete representation of the :ref:`Geometry` and
:ref:`Topology` of a mesh is the key distinguishing characteristic used to
classify a mesh into the different :ref:`MeshTypes`. As a prerequisite to the
proceeding discussion on the taxonomy of the various :ref:`MeshTypes`, this
section provides a high level description of the key constituents of the
:ref:`MeshRepresentation`. Namely, the :ref:`Topology`, :ref:`Geometry` and
:ref:`FieldData` comprising a mesh.

.. _Topology:

Topology
^^^^^^^^

The topology of a mesh, :math:`\mathcal{M}(\Omega) \in \mathbb{R}^d`, is
defined by the collection of topological entities, e.g., the :ref:`Cells`,
:ref:`Faces` and :ref:`Nodes`, comprising the mesh and the associated
*adjacency* information that encodes the topological connections between them,
broadly referred to as :ref:`Connectivity` information.
Each topological entity in the mesh is identified by a unique index, as depicted
in the sample :ref:`UnstructuredMesh` shown in
:numref:`figs/meshRepresentation`. This provides a convenient way to traverse
and refer to individual entities in the mesh.

.. _figs/meshRepresentation:
.. figure:: ../figures/mesh_representation.png
  :align: center
  :scale: 95%
  :alt: Mesh Representation

  Sample unstructured mesh. Each node, cell and face on the mesh has a unique
  index.

In Mint, the three fundamental topological entities comprising a mesh are
(1) :ref:`Cells`, (2) :ref:`Faces`, and (3) :ref:`Nodes`.

.. note::
  The current implementation does not provide first class support for edges and
  associated edge data in 3D. However, this is a feature we are planning
  to support in future versions of Mint.

.. _Cells:

Cells
"""""

A cell, :math:`\mathcal{C}_i`, is given by an ordered list of :ref:`Nodes`,
:math:`\mathcal{C}_i=\{n_0,n_1,...n_k\}`, where each entry,
:math:`n_j \in \mathcal{C}_i`, corresponds to a
unique node index in the mesh. The order of :ref:`Nodes` defining a cell is
determined according to a prescribed local numbering convention [#f2]_ for a
particular cell type. See :numref:`figs/linearCells` and :numref:`figs/q2Cells`.
All Mint :ref:`CellTypes` follow the `CGNS`_ standard local numbering
convention.

.. _Faces:

Faces
"""""

Similarly, a face, :math:`\mathcal{F}_i`, is defined by an ordered list of
:ref:`Nodes`, :math:`\mathcal{F}_i=\{n_0,n_1,...,n_k\}`. Faces are essentially
:ref:`Cells` whose topological dimension is one less than the dimension of the
:ref:`Cells` they are bound to. See :numref:`figs/cellFaces`.
Consequently, the constituent faces of a 3D cell are 2D topological entities,
such as, *triangles* or *quads*, depending on the cell type. The faces of a 2D
cell are 1D topological entities, i.e., *segments*. Last, the faces of a 1D cell
are 0D topological entities, i.e., :ref:`Nodes`.

.. _figs/cellFaces:
.. figure:: ../figures/cell_faces.png
  :align: center
  :alt: Cell Faces

  Constituent faces of a cell in 2D and 3D respectively. the constituent faces
  of a 3D cell are 2D topological entities, such as, *triangles* or *quads*,
  depending on the cell type. The faces of a 2D cell are 1D topological
  entities, i.e., *segments*.

.. admonition:: The Three Face Types

  A mesh face can be bound to either *one* or *two* :ref:`Cells`:

  * :ref:`Faces` bound to two :ref:`Cells`, within the same domain, are called
    **internal faces**.

  * :ref:`Faces` bound to two :ref:`Cells`, across different domains
    (or partitions), are called **internal boundary faces**. Internal boundary
    faces define the communication boundaries where ghost data is exchanged
    between domains.

  * :ref:`Faces` bound to a single cell are called **external boundary faces**.
    External boundary faces (and/or their consistuent nodes) are typically
    associated with a boundary condition.

.. topic:: Face Orientation

   As with :ref:`Cells`, the ordering of the constituent nodes of a face is
   determined by the cell type. However, by convention, the orientation of a
   face is according to an outward pointing face normal, as illustrated in
   :numref:`figs/faceOrientation`.

   .. _figs/faceOrientation:
   .. figure:: ../figures/face_orientation.png
      :align: center
      :alt: Face Orientation

      Face Orientation. (a) From the viewpoint of a cell, its constituent faces
      are oriented according to an outward facing normal. (b) From the viewpoint
      of a face, a face is oriented according to an outward facing normal with
      respect to the first cell abutting to the face, denoted by,
      :math:`\mathcal{C}_0`.

   From the viewpoint of a cell, :math:`\mathcal{C}_k`, its constituent faces,
   defined in the local node numbering of the cell, are oriented
   according to an outward facing normal with respect to the cell,
   :math:`\mathcal{C}_k`. For example, in :numref:`figs/faceOrientation` (a),
   the triangle, :math:`\mathcal{C}_k`, has three faces that are oriented
   according to an outward facing normal and defined using local node numbers
   with respect to their cell as follows, :math:`\mathcal{F}_0=\{0,1\}`,
   :math:`\mathcal{F}_1=\{1,2\}` and :math:`\mathcal{F}_2=\{2,0\}`

   As noted earlier, a face can have at most two adjacent :ref:`Cells`,
   denoted by :math:`\mathcal{C}_0` and :math:`\mathcal{C}_1`. By convention,
   from the viewpoint of a face, :math:`\mathcal{F}_k`, defined using global node
   numbers, the face is oriented according to an outward facing normal with
   respect to the cell corresponding to :math:`\mathcal{C}_0`. As depicted in
   :numref:`figs/faceOrientation` (b), the face denoted by :math:`\mathcal{F}_k`
   has an outward facing normal with respect to :math:`\mathcal{C}_0`,
   consequently it is defined as follows, :math:`\mathcal{F}_k=\{1,2\}`.

.. note::

    By convention, :math:`\mathcal{C}_1` is set to :math:`-1` for
    *external boundary faces*, which are bound to a single cell.

.. _Nodes:

Nodes
"""""

The :ref:`Nodes` are *zero* dimensional topological entities and hence, are the
lowest dimensional constituent entities of a mesh. The :ref:`Nodes` are
associated with the spatial coordinates of the mesh and are used in defining
the topology of the higher dimensional topological entities comprising the mesh,
such as, the :ref:`Cells`, :ref:`Faces`, etc., as discussed earlier.
In a sense, the :ref:`Nodes` provide the means to link the :ref:`Topology`
of the mesh to its constituent :ref:`Geometry` and thereby instantiate the mesh
in physical space.

.. admonition:: Definition

   A mesh node, :math:`\mathcal{n_i}`, is associated with a point,
   :math:`p_i \in \mathbb{R}^d` and provides the means to:

   #. Link the :ref:`Topology` of the mesh to its constituent :ref:`Geometry`

   #. Support one or more *degrees of freedom*, evaluated at the given node
      location.


Notably, the nodes of a mesh are not necessarily just the *vertices* of the mesh.
As discussed in the :ref:`PreliminaryConcepts` section, a mesh is a
discretization of a PDE. Recall, the primary purpose of the mesh is to define
the dicsrete locations, in both *space* and *time*, at which the
*unknown variables* or *degrees of freedom* of the governing PDE are evaluated.
Depending on the numerical scheme employed and the :ref:`CellTypes` used, the
:ref:`Nodes` of a mesh may also be located at cell, face and edge centroids.
For example, in the Finite Element Method (FEM), the nodes for the linear
Lagrange Finite Elements, see :numref:`figs/linearCells`, are located at the
cell *vertices*. However, for quadratic :ref:`CellTypes`,
see :numref:`figs/q2cells`, the *Lagrange* :math:`P^2` finite element,
for the quadrilateral and hexahedron (in 3D) cells, includes as :ref:`Nodes`,
the cell, face and edge (in 3D) centroids in addition to the cell *vertices*.
Other higher order finite elements may involve additional nodes for each edge
and face as well as in the interior of the cell.

.. _Connectivity:

Connectivity
""""""""""""
The topological connections or *adjecencies* between the :ref:`Cells`,
:ref:`Faces` and :ref:`Nodes` comprising the mesh, give rise to a hierarchical
topological structure, depicted in :numref:`figs/topological_structure`, that is
broadly referred to as :ref:`Connectivity` information. At the top level, a
mesh consists of one or more :ref:`Cells`, which constitute the highest
dimensional entity comprising the mesh. Each cell is bounded by zero or more
:ref:`Faces`, each of which is bounded by one or more :ref:`Nodes`.

.. _figs/topological_structure:
.. figure:: ../figures/topological_structure.png
  :align: center
  :scale: 95%
  :alt: Topological Structure

  Hierarchical topological structure illustrating the *downward* and *upward*
  topological connections of the constituent mesh entities supported in Mint.

The topological connections between the constituent entities of the mesh can be
distinguished in (a) *downward* and (b) *upward* topological connections, as
illustrated in :numref:`figs/topological_structure`.

    * The downward topological connections encode the connections from
      higher dimensional mesh entities to lower dimensional entities,
      such as, *cell-to-node*, *face-to-node* or *cell-to-face*.

    * The upward topological connections, also called *reverse connectivities*,
      encode the connections from lower dimensional mesh entities to higher
      dimensional entities, such as, *face-to-cell*.

Two key guiding considerations in the design and implementation of mesh data
structures are *storage* and *computational efficiency*.  In that respect,
the various :ref:`MeshTypes` offer different advantages and tradeoffs. For
example, the inherent regular topology of a :ref:`StructuredMesh` implicitly
defines the :ref:`Connectivity` information. Consequently, the topological
connections between mesh entities can be efficiently computed on-the-fly.
However, for an :ref:`UnstructuredMesh`, the :ref:`Connectivity` information has
to be extracted and stored explicitly so that it is readily available for
computation.

An :ref:`UnstructuredMesh` representation that explicitly stores all :math:`0`
to :math:`d` topological entities and associated *downward* and *upward*
:ref:`Connectivity` information is said to be a *full mesh representation*.
Otherwise, it is called a *reduced mesh representation*. In practice,
it can be prohibitively expensive to store a *full mesh representation*.
Consequently, most applications keep a *reduced mesh representation*.

The question that needs to be addressed at this point is what
:ref:`Connectivity` information is generally required. The answer can vary
depending on the application. The type of operations performed on the mesh
impose the requirements for the :ref:`Connectivity` information needed.
The *minimum sufficient* representation for an :ref:`UnstructuredMesh` is
the *cell-to-node* :ref:`Connectivity`, since, all additional
:ref:`Connectivity` information can be subsequently computed based on
this information.

In an effort to balance both flexibility and simplicity, Mint, in its simplest
form, employs the *minumum sufficient* :ref:`UnstructuredMesh` representation,
consisting of the *cell-to-node* :ref:`Connectivity`. This allows applications
to employ a fairly *light-weight* mesh representation when possible. However,
for applications that demand additional :ref:`Connectivity` information, Mint
provides methods to compute the needed additional information.

.. warning::

   The present implementation of Mint provides first class support for
   *cell-to-node*, *cell-to-face*, *face-to-cell* and *face-to-node*
   :ref:`Connectivity` information for all the :ref:`MeshTypes`.
   Support for additional :ref:`Connectivity` information will be added
   in future versions based on demand by applications.

.. _Geometry:

Geometry
^^^^^^^^
The :ref:`Geometry` of a mesh is given by a set of :ref:`Nodes`.
Let :math:`\mathcal{N}=\{n_0, n_1, n_2, ..., n_k\}` be the finite set of nodes
comprising a mesh, :math:`\mathcal{M}(\Omega) \in \mathbb{R}^d`, where :math:`d`
is the spatial dimension, :math:`d \in \{1,2,3\}`. Each node,
:math:`n_i \in \mathcal{N}`, corresponds to a point,
:math:`p_i \in \mathbb{R}^d`, whose spatial coordinates, i.e., an ordered tuple,
define the physical location of the node in space,
:math:`n_i \in \mathbb{R}^d` . The :ref:`Nodes` link the
:ref:`Geometry` of the mesh to its :ref:`Topology`. The :ref:`Geometry` and
:ref:`Topology` of the mesh collectively define the physical *shape*, *size*
and *location* of the mesh in space.

.. _FieldData:

Field Data
^^^^^^^^^^

The :ref:`FieldData` are used to define various physical quantities over
the constituent mesh entities, i.e., the :ref:`Cells`, :ref:`Faces` and
:ref:`Nodes` of the mesh. Each constituent mesh entity can be associated with
zero or more *fields*, each of which may correspond to a *scalar*, *vector* or
*tensor* quantity, such as, temperature, velocity, pressure, etc.
Essentially, the :ref:`FieldData` are used to define the solution to the
unknown variables of the governing PDE that are evaluated on a given mesh,
as well as, any other auxiliary variables or derived quantities that an
application may need.

.. warning::

   The present implementation of Mint supports :ref:`FieldData` defined on
   :ref:`Cells`, :ref:`Faces` and :ref:`Nodes`. Support for storing
   :ref:`FieldData` on edges will be added in future versions based on
   application demand.


.. #############################################################################
..  MESH TYPES
.. #############################################################################

.. _MeshTypes:

Mesh Types
-----------

The underlying, concrete, representation of the constituent :ref:`Geometry`
and :ref:`Topology` of a mesh is the key defining charachteristic
used in classifying a mesh into the different :ref:`MeshTypes`.
The :ref:`Geometry` and :ref:`Topology` of a mesh is specified in one of the
following three representations:

#. **Implicit Representation**: based on mesh metadata

#. **Explicit Representation**: employs explicitly stored information.

#. **Semi-Implicit Representation**: combines mesh metadata and explicitly
   stored information.

The possible representation combinations of the constituent :ref:`Geometry` and
:ref:`Topology` comprising a mesh define a taxonomy of :ref:`MeshTypes`
summarized in the table below.

.. |structured|   replace:: **Structured Mesh**
.. |curvilinear|  replace:: :ref:`CurvilinearMesh`
.. |rectilinear|  replace:: :ref:`RectilinearMesh`
.. |uniform|      replace:: :ref:`UniformMesh`
.. |unstructured| replace:: :ref:`UnstructuredMesh`
.. |particles|    replace:: :ref:`ParticleMesh`
.. |implicit|     replace:: *implicit*
.. |explicit|     replace:: *explicit*
.. |semi|         replace:: *semi-implicit*

.. raw:: html

      <center>

+------------------+------------+------------+
|   Mesh Type      | Geometry   | Topology   |
|                  |            |            |
+==================+============+============+
| |curvilinear|    | |explicit| | |implicit| |
+------------------+------------+------------+
| |rectilinear|    | |semi|     | |implicit| |
+------------------+------------+------------+
| |uniform|        | |implicit| | |implicit| |
+------------------+------------+------------+
| |unstructured|   | |explicit| | |explicit| |
+------------------+------------+------------+
| |particles|      | |explicit| | |implicit| |
+------------------+------------+------------+

.. raw:: html

      </center>


A brief overview of the distinct characteristics of each of the :ref:`MeshTypes`
is provided in the following sections.

.. #############################################################################
..  STRUCTURED MESH
.. #############################################################################

.. _StructuredMesh:

Structured Mesh
^^^^^^^^^^^^^^^^

A :ref:`StructuredMesh` discretization is characterized by its *ordered*,
*regular*, :ref:`Topology`. A :ref:`StructuredMesh` divides the computational
domain into :ref:`Cells` that are logically arranged on a *regular grid*. The
regular grid topology allows for the constituent :ref:`Nodes`, :ref:`Cells`
and :ref:`Faces` of the mesh to be identified using an *IJK* ordering scheme.

.. admonition:: Numbering and Ordering Conventions in a Structured Mesh

  The *IJK* ordering scheme employs indices along each dimension, typically
  using the letters *i,j,k* for the 1st, 2nd and 3rd dimension respectively.
  The IJK indices can be thought of as counters. Each index counts the number
  of :ref:`Nodes` or :ref:`Cells` along a given dimension. As noted in
  the general :ref:`MeshRepresentation` section, the constituent entities of
  the mesh :ref:`Topology` are associated with a unique index.
  Therefore, a convention needs to be established for mapping the IJK indices to
  the corresponding unique index and *vice-versa*.

  The general convention and what Mint employs is the following:

  * All :ref:`Nodes` and :ref:`Cells` of a :ref:`StructuredMesh` are indexed
    first along the *I*-direction, then along the *J*-direction and last along
    the *K*-direction.

  * Likewise, the :ref:`Faces` of a :ref:`StructuredMesh` are indexed by first
    counting the :ref:`Faces` of each of the :ref:`Cells` along the
    *I*-direction (*I-Faces*), then the *J*-direction (*J-Faces*) and
    last the *K*-direction (*K-Faces*).

One of the important advantages of a :ref:`StructuredMesh` representation is
that the constituent :ref:`Topology` of the mesh is *implicit*. This enables
a convenient way for computing the :ref:`Connectivity` information
automatically without the need to store this information explicitly. For example,
a 2D cell located at :math:`C=(i,j)`, will always have four face neighbors
given by the following indices:

* :math:`N_0=(i-1,j)`,
* :math:`N_1=(i+1,j)`,
* :math:`N_2=(i,j-1)` and
* :math:`N_3=(i,j+1)`

Notably, the neighboring information follows directly from the *IJK* ordering
scheme and therefore does not need to be stored explicitly.

In addition to the convenience of having automatic :ref:`Connectivity`, the
*IJK* ordering of a :ref:`StructuredMesh` offers one other important advantage
over an :ref:`UnstructuredMesh` discretization. The *IJK* ordering
results in coefficient matrices that are *banded*. This enables the use of
specialized algebraic solvers that rely on the *banded* structure of the matrix
that are generally more efficient.

While a :ref:`StructuredMesh` discretization offers several advantages, there
are some notable tradeoffs and considerations. Chief among them,
is the implied restriction imposed by the regular topology of the
:ref:`StructuredMesh`. Namely, the number of :ref:`Nodes` and :ref:`Cells` on
opposite sides must be matching. This requirement makes *local refinement*
less effective, since grid lines need to span across the entire range along
a given dimension. Moreover, meshing of complex geometries, consisting of
sharp features, is complicated and can lead to degenerate :ref:`Cells` that can be
problematic in the computation. These shortcomings are alleviated to an extent
using a *block-structured meshing strategy* and/or *patch-based AMR*, however,
the fundamental limitations still persist.

All :ref:`StructuredMesh` types have implicit :ref:`Topology`. However, depending
on the underlying, concrete representation of the consituent mesh
:ref:`Geometry`, a :ref:`StructuredMesh` is distinguished into three
subtypes:

#. :ref:`CurvilinearMesh`,
#. :ref:`RectilinearMesh`, and,
#. :ref:`UniformMesh`

The key characteristics of each of theses types is discussed in more detail in
the following sections.

.. _CurvilinearMesh:

Curvilinear Mesh
"""""""""""""""""

The :ref:`CurvilinearMesh`, shown in :numref:`figs/curvilinearMeshExample`, is
logically a *regular* mesh, however, in contrast to the :ref:`RectilinearMesh`
and :ref:`UniformMesh`, the :ref:`Nodes` of a :ref:`CurvilinearMesh` are not
placed along the *Cartesian* grid lines. Instead, the equations of the governing
PDE are transformed from the *Cartesian* coordinates to a new coordinate system,
called a *curvilinear coordinate system*. Consequently, the :ref:`Topology` of
a :ref:`CurvilinearMesh` is *implicit*, however, its :ref:`Geometry`, given
by the constituent :ref:`Nodes` of the mesh, is *explicit*.

.. _figs/curvilinearMeshExample:
.. figure:: ../figures/structured_curvilinear_mesh.png
  :align: center
  :scale: 55%
  :alt: Sample Curvilinear Mesh

  Sample Curvilinear Mesh example.

The mapping of coordinates to the *curvilinear coordinate system*
facilitates the use of structured meshes for bodies of arbitrary shape. Note,
the axes defining the *curvilinear coordinate system* do not need to be straight
lines. They can be curves and align with the contours of a solid body. For this
reason, the resulting :ref:`CurvilinearMesh` is often called a *mapped mesh* or
*body-fitted mesh*.

.. _RectilinearMesh:

Rectilinear Mesh
"""""""""""""""""

A :ref:`RectilinearMesh`, depicted in :numref:`figs/rectilinearMeshExample`,
divides the computational domain into a set of rectangular :ref:`Cells`,
arranged on a *regular lattice*. However, in contrast to the
:ref:`CurvilinearMesh`, the :ref:`Geometry` of the mesh is not mapped to a
different coordinate system. Instead, the rows and columns of :ref:`Nodes`
comprising a :ref:`RectilinearMesh` are parallel to the axis of the *Cartesian*
coordinate system. Due to this restriction, the geometric domain and resulting
mesh are always rectangular.

.. _figs/rectilinearMeshExample:
.. figure:: ../figures/structured_rectilinear_mesh.png
  :align: center
  :scale: 35%
  :alt: Sample Rectilinear Mesh

  Sample Rectilinear Mesh example.

The :ref:`Topology` of a :ref:`RectilinearMesh` is *implicit*, however, its
constituent :ref:`Geometry` is *semi-implicit*. Although, the
:ref:`Nodes` are aligned with the *Cartesian* coordinate axis, the spacing
between adjacent :ref:`Nodes` can vary. This allows a :ref:`RectilinearMesh`
to have tighter spacing over regions of interest and be sufficiently coarse in
other parts of the domain. Consequently, the spatial coordinates of the
:ref:`Nodes` along each axis are specified explicitly in a seperate array
for each coordinate axis, i.e., :math:`x`, :math:`y` and :math:`z` arrays for
each dimension respectively. Given the *IJK* index of a node, its corresponding
physical coordinates can be obtained by taking the *Cartesian* product of the
corresponding coordinate along each coordinate axis. For this reason, the
:ref:`RectilinearMesh` is sometimes called a *product* mesh.

.. _UniformMesh:

Uniform Mesh
"""""""""""""

A :ref:`UniformMesh`, depicted in :numref:`figs/uniformMeshExample`, is the
simplest of all three :ref:`StructuredMesh` types, but also, relatively the most
restrictive of all :ref:`MeshTypes`. As with the :ref:`RectilinearMesh`,
a :ref:`UniformMesh` divides the computational domain into a set of rectangular
:ref:`Cells` arranged on a *regular lattice*. However, a :ref:`UniformMesh`
imposes the additional restriction that :ref:`Nodes` are uniformly distributed
parallel to each axis. Therefore, in contrast to the :ref:`RectilinearMesh`, the
spacing between adjacent :ref:`Nodes` in a :ref:`UniformMesh` is constant.

.. _figs/uniformMeshExample:
.. figure:: ../figures/structured_uniform_mesh.png
  :align: center
  :scale: 35%
  :alt: Sample Uniform Mesh

  Sample Uniform Mesh example.

The inherent constraints of a :ref:`UniformMesh` allow for a more compact
representation. Notably, both the :ref:`Topology` and :ref:`Geometry` of a
:ref:`UniformMesh` are *implicit*. Given the origin of the mesh,
:math:`X_0=(x_0,y_0,z_0)^T`, i.e., the coordinates of the lowest corner of the
rectangular domain, and spacing along each direction, :math:`H=(h_x,h_y,h_z)^T`,
the spatial coordinates of any point, :math:`\hat{p}=(p_x,p_y,p_z)^T`,
corresponding to a node with lattice coordinates, :math:`(i,j,k)`, are
computed as follows:

.. math::
    :nowrap:

    \begin{eqnarray}
      p_x &=& x_0 &+& i &\times& h_x \\
      p_y &=& y_0 &+& j &\times& h_y \\
      p_z &=& z_0 &+& k &\times& h_z \\
    \end{eqnarray}


.. #############################################################################
..  UNSTRUCTURED MESH
.. #############################################################################

.. _UnstructuredMesh:

Unstructured Mesh
^^^^^^^^^^^^^^^^^^

The impetus for an :ref:`UnstructuredMesh` discretization is largely prompted
by the need to model physical phenomena on complex geometries. In relation to
the various :ref:`MeshTypes`, an :ref:`UnstructuredMesh` discretization provides
the most flexibility. Notably, an :ref:`UnstructuredMesh` can accomodate
different :ref:`CellTypes` and does not enforce any constraints or particular
ordering on the constituent :ref:`Nodes` and :ref:`Cells`. This makes an
:ref:`UnstructuredMesh` discretization particularly attractive, especially for
applications that require *local adaptive mesh refinement* [#f3]_ and deal with
complex geometries.

Generally, the advantages of using an :ref:`UnstructuredMesh` come at the cost
of an increase in memory requirements and computational intensity. This is
due to the inherently *explicit*, :ref:`MeshRepresentation`
required for an :ref:`UnstructuredMesh`. Notably, both :ref:`Topology` and
:ref:`Geometry` are represented explicitly thereby increasing the storage
requirements and computational time needed per operation. For example, consider
a stencil operation. For a :ref:`StructuredMesh`, the neighbor indices needed
by the stencil can be automatically computed directly from the *IJK* ordering,
a relatively fast and local operation. However, to obtain the neighbor indices
in an :ref:`UnstructuredMesh`, the arrays that store the associated
:ref:`Connectivity` information need to be accessed, resulting in additional
load/store operations that are generaly slower.

Depending on the application, the constituent :ref:`Topology` of an
:ref:`UnstructuredMesh` may employ a:

#. :ref:`SingleCellTopology`, i.e., consisting of :ref:`Cells` of the *same type*, or,
#. :ref:`MixedCellTopology`, i.e., consisting of :ref:`Cells` of different type, i.e., *mixed cell type*.

There are subtle differrences in the underlying :ref:`MeshRepresentation` that
can result in a more compact and efficient representation when the
:ref:`UnstructuredMesh` employs a :ref:`SingleCellTopology`. The following
sections discuss briefly these differences and other key aspects of the
:ref:`SingleCellTopology` and :ref:`MixedCellTopology` representations.
Moreover, the list of natively supported :ref:`CellTypes` that can be used with
an :ref:`UnstructuredMesh` is presented, as well as, the steps necessary to
:ref:`AddACellType` in Mint.

.. note::
    In an effort to balance both flexibility and simplicity, Mint, in its simplest
    form, employs the *minumum sufficient* :ref:`UnstructuredMesh`
    :ref:`MeshRepresentation`, consisting of the *cell-to-node*
    :ref:`Connectivity`. This allows applications to employ a fairly
    *light-weight* mesh representation when possible. However, for applications
    that demand additional :ref:`Connectivity` information, Mint provides
    methods to compute the needed additional information.


.. _SingleCellTopology:

Single Cell Type Topology
""""""""""""""""""""""""""

An :ref:`UnstructuredMesh` with :ref:`SingleCellTopology` consists of a
collection of :ref:`Cells` of the same cell type. Any :ref:`StructuredMesh`
can be treated as an :ref:`UnstructuredMesh` with :ref:`SingleCellTopology`,
in which case, the resulting :ref:`Cells` would either be *segments* (in 1D),
*quadrilaterals* (in 2D) or *hexahedrons* (in 3D). However, an
:ref:`UnstructuredMesh` can have arbitrary :ref:`Connectivity` and does not
impose any ordering constraints. Moreover, the :ref:`Cells` can also be
*triangular* (in 2D) or *tetrahedral* (in 3D). The choice of cell type generally
depends on the application, the physics being modeled, and the numerical
scheme employed. An example tetrahedral :ref:`UnstructuredMesh` of the F-17
blended wing fuselage configuration is shown in
:numref:`figs/UnstructuredMeshSingleShape`. For this type of complex geometries
it is nearly impossible to obtain a :ref:`StructuredMesh` that is adequate for
computation.

.. _figs/unstructuredMeshSingleShape:
.. figure:: ../figures/f17.png
   :align: center
   :scale: 35%
   :alt: Sample Unstructured Mesh (single shape topology)

   Sample unstructured tetrahedral mesh of the F-17 blended wing fuselage
   configuration.

Mint's :ref:`MeshRepresentation` of an :ref:`UnstructuredMesh` with
:ref:`SingleCellTopology` consists of a the cell type specification and the
cell-to-node :ref:`Connectivity` information. The :ref:`Connectivity` information
is specified with a flat array consisting of the node indices that comprise
each cell. Since the constituent mesh :ref:`Cells` are of the same type,
cell-to-node information for a particular cell can be obtained by accessing
the :ref:`Connectivity` array with a constant stride, where the stride
corresponds to the number of :ref:`Nodes` of the cell type being used. This is
equivalent to a 2D row-major array layout where the number of rows corresponds
to the number of :ref:`Cells` in the mesh and the number of columns corresponds
to the *stride*, i.e., the number of :ref:`Nodes` per cell.

.. _figs/singleCellTypeRep:
.. figure:: ../figures/SingleCellTypeMesh.png
    :align: center
    :alt: Mesh Representation of the Unstructured Mesh with Single Cell Topology

    :ref:`MeshRepresentation` of an :ref:`UnstructuredMesh` with
    :ref:`SingleCellTopology` consiting of *triangular* :ref:`Cells`. Knowing
    the cell type enables traversing the cell-to-node :ref:`Connectivity` array
    with a constant stride of :math:`3`, which corresponds to the number of
    constituent :ref:`Nodes` of each triangle.

This simple concept is best illustrated with an example.
:numref:`figs/singleCellTypeRep` depicts a sample :ref:`UnstructuredMesh` with
:ref:`SingleCellTopology` consisting of :math:`N_c=4` triangular :ref:`Cells`.
Each triangular cell, :math:`C_i`, is defined by :math:`||C_i||` :ref:`Nodes`.
In this case, :math:`||C_i||=3`.

.. note::

    The number of :ref:`Nodes` of the cell type used to define an
    :ref:`UnstructuredMesh` with :ref:`SingleCellTopology`,
    denoted by :math:`||C_i||`, corresponds to the constant stride used to
    access the flat cell-to-node :ref:`Connectivity` array.

Consequently, the length of the cell-to-node :ref:`Connectivity` array
is then given by :math:`N_c \times ||C_i||`. The node indices for each of the
cells are stored from left to right. The base offset for a given cell is given
as a multiple of the cell index and the *stride*. As illustrated in
:numref:`figs/singleCellTypeRep`, the base offset for cell :math:`C_0` is
:math:`0 \times 3 = 0`, the offest for cell :math:`C_1` is
:math:`1 \times 3 = 1`, the offset for cell :math:`C_2` is
:math:`2 \times 3 = 6` and so on.

.. admonition:: Direct Stride Cell Access in a Single Cell Type Topology UnstructuredMesh

    In general, the :ref:`Nodes` of a cell, :math:`C_i`, of an :ref:`UnstructuredMesh`
    with :ref:`SingleCellTopology` and cell stride :math:`||C_i||=k`, can be
    obtained from a given cell-to-node :ref:`Connectivity` array as follows:

    .. math::
      :nowrap:

      \begin{eqnarray}
        n_0 &=& cell\_to\_node[ i \times k     ] \\
        n_1 &=& cell\_to\_node[ i \times k + 1 ] \\
        ... \\
        n_k &=& cell\_to\_node[ i \times k + (k-1)]
      \end{eqnarray}

.. raw:: html

      <center>

+------------------+------------+------------+------------+
|   Cell Type      | Stride     | Topological| Spatial    |
|                  |            | Dimension  | Dimension  |
+==================+============+============+============+
| *Quadrilateral*  |     4      |  2         | 2,3        |
+------------------+------------+------------+------------+
| *Triangle*       |     3      |  2         | 2,3        |
+------------------+------------+------------+------------+
| *Hexahdron*      |     8      |  3         |  3         |
+------------------+------------+------------+------------+
| *Tetrahedron*    |     4      |  3         |  3         |
+------------------+------------+------------+------------+

.. raw:: html

      </center>

The same procedure follows for any cell type. Thereby, the stride for a mesh
consisting of *quadrilaterals* is :math:`4`, the stride for a mesh consisting
of *tetrahedrons* is :math:`4` and the stride for a mesh consisting of
*hexahedrons* is :math:`8`. The table above summarizes the possible
:ref:`CellTypes` that can be employed for an :ref:`UnstructuredMesh` with
:ref:`SingleCellTopology`, corresponding *stride* and applicalble topological
and spatial dimension.


.. _MixedCellTopology:

Mixed Cell Type Topology
"""""""""""""""""""""""""

An :ref:`UnstructuredMesh` with :ref:`MixedCellTopology` provides the most
flexibility relative to the other :ref:`MeshTypes`. Similar to the
:ref:`SingleCellTopology` :ref:`UnstructuredMesh`,  the constituent
:ref:`Nodes` and :ref:`Cells` of a :ref:`MixedCellTopology` :ref:`UnstructuredMesh`
can have arbitrary ordering. Both :ref:`Topology` and :ref:`Geometry` are
*explicit*. However, a :ref:`MixedCellTopology` :ref:`UnstructuredMesh` may
consist :ref:`Cells` of different cell type. Hence, the cell topology and
cell type is said to be *mixed*.

.. note::
   The constituent :ref:`Cells` of an :ref:`UnstructuredMesh` with
   :ref:`MixedCellTopology` have a *mixed cell type*. For this reason,
   an :ref:`UnstructuredMesh` with :ref:`MixedCellTopology` is sometimes also
   called a *mixed cell mesh* or *hybrid mesh*.

.. _figs/unstructuredMeshMixedShape:
.. figure:: ../figures/unstructured_mixed_mesh.png
   :align: center
   :scale: 95%
   :alt: Sample Unstrucrured Mesh (mixed shape topology)

   Sample :ref:`UnstructuredMesh` with :ref:`MixedCellTopology` of a Generic
   wing/fuselage configuration. The mesh consists of high-aspect ratio prism
   cells in the viscous region of the computational domain to accurately capture
   the high gradients across the boundary layer and tetrahedra cells for the
   inviscid/Euler portion of the mesh.

Several factors need to be taken in to account when selecting
the cell topology of the mesh. The physics being modeled, the PDE
discretization employed and the required simulation fidelity are chief among
them. Different :ref:`CellTypes` can have superior properties for certain
calculations. The continuous demand for increasing fidelity in physics-based
predictive modeling applications has prompted practitioners to employ a
:ref:`MixedCellTopology` :ref:`UnstructuredMesh` discretization in order to
accurately capture the underlying physical phenomena.

For example, for Navier-Stokes *viscous* fluid-flow computations, at high
Reynolds numbers, it is imperative to capture the high gradients across the
boundary layer normal to the wall. Typically, high-aspect ratio, anisotropic
*triangular prisms* or *hexahedron* :ref:`Cells` are used for discretizing the
viscous region of the computational domain, while isotropic *tetrahedron* or
*hexahedron* :ref:`Cells` are used in the *inviscid* region to solve the Euler
equations. The sample :ref:`MixedCellTopology` :ref:`UnstructuredMesh`, of a
Generic Wing/Fuselage configuration, depicted in
:numref:`figs/unstructuredMeshMixedShape`, consists of *triangular prism*
:ref:`Cells` for the *viscous* boundary layer portion of the domain that are
stitched to *tetrahedra* :ref:`Cells` for the inviscid/Euler portion of the
mesh.

The added flexibility enabled by employing a :ref:`MixedCellTopology`
:ref:`UnstructuredMesh` imposes additional requirements to the underlying
:ref:`MeshRepresentation`. Most notably, compared to the
:ref:`SingleCellTopology` :ref:`MeshRepresentation`, the cell-to-node
:ref:`Connectivity` array can consist :ref:`Cells` of different cell type,
where each cell can have a different number of :ref:`Nodes`. Consequently, the
simple stride array access indexing scheme, used for
the :ref:`SingleCellTopology` :ref:`MeshRepresentation`, cannot be employed to
obtain cell-to-node information. For a :ref:`MixedCellTopology` an
*indirect addressing* access scheme must be used instead.

.. _figs/mixedCellTypeRep:
.. figure:: ../figures/MixedCellTypeMesh.png
    :align: center
    :alt: Mesh Representation of the Unstructured Mesh with Mixed Cell Topology

    :ref:`MeshRepresentation` of a :ref:`MixedCellTopology`
    :ref:`UnstructuredMesh` with a total of :math:`N=3` :ref:`Cells`,
    :math:`2` *triangles* and :math:`1` *quadrilateral*. The
    :ref:`MixedCellTopology` representation consists of two additional arrays.
    First, the *Cell Offsets* array, an array of size :math:`N+1`, where the
    first :math:`N` entries store the starting position to the flat
    cell-to-node :ref:`Connectivity` array for each cell. The last entry of
    the *Cell Offsets* array stores the total length of the :ref:`Connectivity`
    array. Second, the :ref:`CellTypes` array , an array of size :math:`N`,
    which stores the cell type of each constituent cell of the mesh.


There are a number of ways to represent a :ref:`MixedCellTopology` mesh.
In addition to the cell-to-node :ref:`Connectivity` array, Mint's
:ref:`MeshRepresentation` for a :ref:`MixedCellTopology` :ref:`UnstructuredMesh`
employs two additional arrays. See sample mesh and corresponding
:ref:`MeshRepresentation` in :numref:`figs/mixedCellTypeRep`.
First, the *Cell Offsets* array is used to provide indirect addressing to
the cell-to-node information of each constituent mesh cell. Second, the
:ref:`CellTypes` array is used to store the cell type of each cell in the
mesh.

The *Cell Offsets* is an array of size :math:`N+1`, where the first
:math:`N` entries, corresponding to each cell in the mesh, store the
start index position to the cell-to-node :ref:`Connectivity` array for the
given cell. The last entry of the *Cell Offsets* array stores the total
length of the :ref:`Connectivity` array. Moreover, the number of constituent
cell :ref:`Nodes` for a given cell can be directly computed by subtracting a
Cell's start index from the next adjacent entry in the  *Cell Offsets* array.

However, knowing the start index position to the cell-to-node
:ref:`Connectivity` array and number of constituent :ref:`Nodes` for a given
cell is not sufficient to disambiguate and deduce the cell type. For example,
both *tetrahedron* and *quadrilateral* :ref:`Cells` are defined by :math:`4`
:ref:`Nodes`. The cell type is needed in order to correctly interpret the
:ref:`Topology` of the cell according to the cell's local numbering.
Consequently, the :ref:`CellTypes` array, whose length is :math:`N`,
corresponding to the number of cells in the mesh, is used to store the cell
type for each constituent mesh cell.

.. admonition:: Indirect Address Cell Access in a Mixed Cell Type Topology UnstructuredMesh

    In general, for a given cell, :math:`C_i`, of a :ref:`MixedCellTopology`
    :ref:`UnstructuredMesh`, the number of :ref:`Nodes` that define the
    cell, :math:`||C_i||`, is given by:

    .. math::
      :nowrap:

      \begin{eqnarray}
        k = ||C_i|| &=& cells\_offset[ i + 1 ] - cells\_offset[ i ] \\
      \end{eqnarray}

    The corresponding cell type is directly obtained from the :ref:`CellTypes`
    array:

    .. math::
      :nowrap:

      \begin{eqnarray}
        ctype &=& cell\_types[ i ] \\
      \end{eqnarray}

    The list of constituent cell :ref:`Nodes` can then obtained from the
    cell-to-node :ref:`Connectivity` array as follows:

    .. math::
      :nowrap:

      \begin{eqnarray}

        offset &=& cells\_offset[ i+1 ] \\
        k      &=& cells\_offset[ i + 1 ] - cell\_offset[ i ] \\

        \\

        n_0 &=& cell\_to\_node[ offset     ] \\
        n_1 &=& cell\_to\_node[ offset + 1 ] \\
        ... \\
        n_k &=& cell\_to\_node[ offset + (k-1)]
      \end{eqnarray}


.. _CellTypes:

Cell Types
"""""""""""

Mint currently supports the common Linear :ref:`CellTypes`,
depicted in :numref:`figs/linearCells`, as well as, support for high-order,
quadratic, quadrilateral and hexahedron :ref:`Cells`, see :numref:`figs/q2Cells`.

.. _figs/linearCells:
.. figure:: ../figures/linear_cell_types.png
  :align: center
  :scale: 95%
  :alt: Supported linear cell types.

  List of supported linear cell types and their respective local node
  numbering.

.. _figs/q2Cells:
.. figure:: ../figures/q2_cell_types.png
  :align: center
  :scale: 95%
  :alt: Supported high-order cell types.

  List of supported quadratic cell types and their respective local node
  numbering.

.. note::

  All Mint :ref:`CellTypes` follow the `CGNS`_ standard local node numbering
  conventions.

Moreover, Mint is designed to be extensible. It is relatively straightforward
to :ref:`AddACellType` in Mint. Each of the :ref:`CellTypes` in Mint simply
encode the following attributes:

* the cell's topology, e.g., number of nodes, faces, local node numbering etc.,
* the corresponding VTK type, used for VTK dumps, and,
* the associated blueprint name, conforming to the `Blueprint`_
  conventions, used for storing the mesh in `Sidre`_

.. warning::
   The `Blueprint`_ specification does not currently support the following
   cell types:

   #. Transitional cell types,  Pyramid(``mint::PYRAMID``)
      and Prism(``mint::PRISM``)

   #. Quadratic cells, the 9-node, quadratic Quadrilateral(``mint::QUAD9``) and
      the 27-node, quadratic Hexahedron(``mint::HEX27``)

.. _AddACellType:

Add a New Cell Type
"""""""""""""""""""

.. warning::
   This section is under construction.

.. #############################################################################
..  PARTICLE MESH
.. #############################################################################

.. _ParticleMesh:

Particle Mesh
^^^^^^^^^^^^^^

A :ref:`ParticleMesh`, depicted in :numref:`figs/particleMesh`, discretizes the
computational domain by a set of *particles* which correspond to the :ref:`Nodes`
at which the solution is evaluated. A :ref:`ParticleMesh` is commonly employed in
the so called *particle* methods, such as, *Smoothed Particle Hydrodynamics*
(SPH) and *Particle-In-Cell* (PIC) methods, which are used in a variety
of applications ranging from astrophysics and cosmology simulations to plasma
physics.

There is no special ordering imposed on the particles. Therefore, the particle
coordinates are explicitly specified by nodal coordinates, similar to an
:ref:`UnstructuredMesh`. However, the particles are not connected to form a
*control volume*, i.e., a filled region of space. Consequently,
a :ref:`ParticleMesh` does not have :ref:`Faces` and any associated
:ref:`Connectivity` information.  For this reason, methods that
employ a :ref:`ParticleMesh` discretization are often referred to as
*meshless* or *mesh-free* methods.

.. _figs/particleMesh:
.. figure:: ../figures/particles.png
   :align: center
   :scale: 35%
   :alt: Sample Particle Mesh

   Sample :ref:`ParticleMesh` within a box domain.

A :ref:`ParticleMesh` can be thought of as having *explicit* :ref:`Geometry`,
but, *implicit* :ref:`Topology`. Mint's :ref:`MeshRepresentation` for a
:ref:`ParticleMesh`, associates the constituent particles with the :ref:`Nodes`
of the mesh. The :ref:`Nodes` of the :ref:`ParticleMesh` can
also be thought of as :ref:`Cells` that are defined by a single node index.
However, since this information can be trivially obtained there is no need to
be stored explicitly.

.. note::

    A :ref:`ParticleMesh` can only store variables at its constituent
    particles, i.e., the :ref:`Nodes` of the mesh. Consequently,
    a :ref:`ParticleMesh` in Mint can only be associated with node-centered
    :ref:`FieldData`.

.. #############################################################################
..  COMPONENT ARCHITECTURE
.. #############################################################################

.. _Architecture:

Overarching Component Architecture
-----------------------------------

This section links the core concepts, presented in the
:ref:`MeshRepresentation` and :ref:`MeshTypes` sections, to the underlying
implementation of the Mint :ref:`MeshDataModel`.
The :ref:`Architecture` of Mint's :ref:`MeshDataModel` consists of a class
hierarchy that follows directly the taxonomy of :ref:`MeshTypes` discussed
earlier. The constituent classes of the :ref:`MeshDataModel` are combined
using a mix of class *inheritance* and *composition*, as illustrated in
the class diagram depicted in :numref:`figs/classDiagram`.

.. _figs/classDiagram:
.. figure:: ../figures/class_diagram.png
   :align: center
   :scale: 50%
   :alt: Mint Class Hierarchy Diagram

   :ref:`Architecture` of the Mint :ref:`MeshDataModel`, depicting the
   core mesh classes and the inter-relationship between them.
   The solid arrows indicate an *inheritance* relationship, while the dashed
   arrows indicate an *ownership* relationship between two classes.

At the top level, :ref:`TheMeshBaseClass`, implemented in ``mint::Mesh``,
stores common mesh attributes and fields. Moreover, it defines a unified
Application Programming Interface (API) for the various :ref:`MeshTypes`. See the
`Mint Doxygen API Documentation`_ for a complete specification of
the API. The :ref:`ConcreteMeshClasses` extend :ref:`TheMeshBaseClass`
and implement the :ref:`MeshRepresentation` for each of the
:ref:`MeshTypes` respectively. The ``mint::ConnectivityArray`` and
``mint::MeshCoordinates`` classes, are the two main internal support classes
that underpin the implementation of the :ref:`ConcreteMeshClasses` and facilitate
the representation of the constituent :ref:`Geometry` and :ref:`Topology`
of the mesh.

.. note::

    All Mint classes and functions are encapsulated in the ``axom::mint``
    namespace.

.. _TheMeshBaseClass:

The Mesh Base Class
^^^^^^^^^^^^^^^^^^^^

:ref:`TheMeshBaseClass` stores common attributes associated with a mesh.
Irrespective of the mesh type, a Mint mesh has two identifiers. The mesh
*BlockID* and mesh *DomainID*, which are assigned by domain decomposition.
Notably, the computational domain can consist of one or more blocks, which are
usually defined by the user or application. Each block is then subsequently
partitioned to multiple domains that are distributed across processing units
for parallel computation. For example, a sample block and domain decomposition
is depicted in :numref:`figs/decomp`. Each of the constituent domains is
represented by a corresponding ``mint::Mesh`` instance, which in aggregate
define the entire problem domain.

.. _figs/decomp:
.. figure:: ../figures/decomp.png
   :align: center
   :scale: 35%
   :alt: Block and Domain Decomposition.

   Sample block & domain decomposition of the computational domain.
   The computational domain is defined using :math:`3` blocks (left). Each block
   is further partitioned into two or more domains(right). A ``mint::Mesh``
   instance represents one of the constituent domains used to define the
   overall problem domain.

.. note::

  A ``mint::Mesh`` instance provides the means to store the mesh
  *BlockID* and *DomainID* respectively. However, Mint does not impose a
  numbering or partitioning scheme. Assignment of the *BlockID* and *DomainID*
  is handled at the application level and by the underlying mesh partitioner
  that is being employed.

Moreover, each ``mint::Mesh`` instance has associated :ref:`MeshFieldData`,
represented by  the ``mint::FieldData`` class. Each of the constituent
topological mesh entities, i.e., the :ref:`Cells`, :ref:`Faces` and :ref:`Nodes`
comprising the mesh, has a handle to a corresponding ``mint::FieldData``
instance. The ``mint::FieldData`` object essentialy provides a container to
store and manage a collection of fields, defined over the corresponding mesh
entity.

.. warning::

   Since a :ref:`ParticleMesh` is defined by a set of :ref:`Nodes`, it
   can only store :ref:`FieldData` at its constituent :ref:`Nodes`. All
   other supported :ref:`MeshTypes` can have :ref:`FieldData` associated with
   their constituent :ref:`Cells`, :ref:`Faces` and :ref:`Nodes`.

.. _MeshFieldData:

Mesh Field Data
^^^^^^^^^^^^^^^^^

A ``mint::FieldData`` instance typically stores multiple fields.
Each field is represented by an instance of a ``mint::Field`` object
and defines a named numerical quantity, such as, *mass*, *velocity*,
*temperature*, etc., defined on a given mesh. Moreover, a field can be either
*single-component*, i.e., a *scalar* quantity, or, *multi-component*, e.g.,
a *vector* or *tensor* quantity. Typically, a field represents some physical
quantity that is being modeled, or, an auxiliary quantity that is needed to
perform a particular calculation.

In addition, each ``mint::Field`` instance can be of different data type.
The ``mint::FieldData`` object can store different types of fields.
For example, *floating point* quantities i.e., ``float`` or ``double``,
as well as, *integral* quantities, i.e., ``int32_t``, ``int64_t``, etc. This is
accomplished using a combination of C++ templates and *inheritance*. The
``mint::Field`` object is an abstract base class that defines a type-agnostic
interface to encapsulate a field. Since ``mint::Field`` is an abstract
base class, it is not instantiated directly. Instead, all fields are created by
instantiating a ``mint::FieldVariable`` object, a class templated on data type,
that derives from the ``mint::Field`` base class. For example, the code snippet
below illustrates how fields of different type can be instantiated.

.. code-block:: cpp

    ...

    // create a scalar field to store mass as a single precision quantity
    mint::Field* mass = new mint::FieldVariable< float >( "mass", size );

    // create a velocity vector field as a double precision floating point quantity
    constexpr int NUM_COMPONENTS = 3;
    mint::Field* vel = new mint::FieldVariable< double >( "vel", size, NUM_COMPONENTS );

    ...

Generally, in application code, it is not necessary to create fields using the
``mint::FieldVariable`` class directly. The ``mint::Mesh`` object provides
convenience methods for adding, removing and accessing fields on a mesh.
Consult the :ref:`sections/tutorial` for more details on
:ref:`workingWithFieldDataOnAMesh`.

.. _ConcreteMeshClasses:

Concrete Mesh Classes
^^^^^^^^^^^^^^^^^^^^^^

The :ref:`ConcreteMeshClasses`, extend :ref:`TheMeshBaseClass` and implement
the underlying :ref:`MeshRepresentation` of the various :ref:`MeshTypes`,
depicted in :numref:`figs/meshtypes`.

.. _figs/meshtypes:
.. figure:: ../figures/meshtypes.png
   :align: center
   :scale: 35%
   :alt: Supported Mesh Types.

   Depiction of the supported :ref:`MeshTypes` with labels of the corresponding
   Mint class used for the underlying :ref:`MeshRepresentation`.


Structured Mesh
""""""""""""""""

All :ref:`StructuredMesh` types in Mint can be represented by an instance of
the ``mint::StructuredMesh`` class, which derives directly from :ref:`TheMeshBaseClass`,
``mint::Mesh``. The ``mint::StructuredMesh`` class is also an abstract base class
that encapsulates the implementation of the *implicit*, *ordered* and *regular*
:ref:`Topology` that is common to all :ref:`StructuredMesh` types. The
distinguishing characteristic of the different :ref:`StructuredMesh` types is
the representation of the constituent :ref:`Geometry`. Mint implements each of
the different :ref:`StructuredMesh` types by a corresponding class, which
derives directly from ``mint::StructuredMesh`` and thereby inherit its
implicit :ref:`Topology` representation.

Consequently, support for the :ref:`UniformMesh` is implemented in
``mint::UniformMesh``. The :ref:`Geometry` of a :ref:`UniformMesh` is
*implicit*, given by two attributes, the mesh *origin* and *spacing*.
Consequently, the ``mint::UniformMesh`` consists of two data members to
store the *origin* and *spacing* of  the :ref:`UniformMesh` and provides
functionality for evaluating the spatial coordinates of a node given its
corresponding IJK lattice coordinates.

The following code snippet provides a simple example illustrating how to
construct and operate on a 2-D :ref:`UniformMesh`.

.. literalinclude:: ../../../examples/mint_uniform_mesh.cpp
   :language: C++
   :linenos:

Similarly, support for the :ref:`RectilinearMesh` is implemented in
``mint::RectilinearMesh``. The constituent :ref:`Geometry` representation of
the :ref:`RectilinearMesh` is *semi-implicit*. The spatial coordinates
of the :ref:`Nodes` along each axis are specified explicitly while the
coordinates of the interior :ref:`Nodes` are evaluated by taking the
*Cartesian* product of the corresponding coordinate along each coordinate
axis. The ``mint::RectilinearMesh`` consists of seperate arrays to
store the coordinates along each axis for the *semi-implicit* :ref:`Geometry`
representation of the :ref:`RectilinearMesh`.

The following code snippet provides a simple example illustrating how to
construct and operate on a 2-D :ref:`RectilinearMesh`.

.. literalinclude:: ../../../examples/mint_rectilinear_mesh.cpp
   :language: C++
   :linenos:

Support for the :ref:`CurvilinearMesh` is implemented by the
``mint::CurvilinearMesh`` class. The :ref:`CurvilinearMesh` requires *explicit*
representation of its constituent :ref:`Geometry`. The ``mint::CurvilinearMesh``
makes use of the ``mint::MeshCoordinates`` class to explicitly represent the
spatial coordinates associated with the constituent :ref:`Nodes` of the mesh.

The following code snippet provides a simple example illustrating how to
construct and operate on a 2-D :ref:`CurvilinearMesh`.

.. literalinclude:: ../../../examples/mint_curvilinear_mesh.cpp
   :language: C++
   :linenos:

Unstructured Mesh
""""""""""""""""""

Mint's :ref:`UnstructuredMesh` representation is provided by the
``mint::UnstructuredMesh`` class, which derives directly from the
:ref:`TheMeshBaseClass`, ``mint::Mesh``. An :ref:`UnstructuredMesh` has both
*explicit* :ref:`Geometry` and :ref:`Topology`. As with the
``mint::CurvilinearMesh`` class, the *explicit* :ref:`Geometry` representation
of the :ref:`UnstructuredMesh` employs the ``mint::MeshCoordinates``. The
constituent :ref:`Topology` is handled by the ``mint::ConnectivityArray``,
which is employed for the representation of all the topological
:ref:`Connectivity` information, i.e., *cell-to-node*, *face-to-node*,
*face-to-cell*, etc.

.. note::

  Upon construction, a ``mint::UnstructuredMesh`` instance consists of the
  *minimum sufficient* representation for an :ref:`UnstructuredMesh` comprised
  of the *cell-to-node* :ref:`Connectivity` information.
  Applications that require face :ref:`Connectivity` information must
  explicitly call the ``initializeFaceConnectivity()`` method on the
  corresponding :ref:`UnstructuredMesh` object.


Depending on the cell :ref:`Topology` being employed, an :ref:`UnstructuredMesh`
can be classified as either a :ref:`SingleCellTopology` :ref:`UnstructuredMesh`
or a :ref:`MixedCellTopology` :ref:`UnstructuredMesh`. To accomodate these
two different representations, the ``mint::UnstructuredMesh`` class, is templated
on ``CELL_TOPOLOGY``. Internally, the template argument is used to indicate
the type of ``mint::ConnectivityArray`` to use, i.e., whether,
*stride access addressing* or *indirect addressing* is used, for
:ref:`SingleCellTopology` and :ref:`MixedCellTopology` respectively.


Particle Mesh
""""""""""""""

Support for the :ref:`ParticleMesh` representation is implemented in
``mint::ParticleMesh``, which derives directly from :ref:`TheMeshBaseClass`,
``mint::Mesh``. A :ref:`ParticleMesh` discretizes the domain by a set
of particles, which correspond to the constituent :ref:`Nodes` of the mesh.
The :ref:`Nodes` of a :ref:`ParticleMesh` can also be thought of as :ref:`Cells`,
however, since this information is trivially obtrained, there is not need
to be stored explicitly, e.g., using a :ref:`SingleCellTopology`
:ref:`UnstructuredMesh` representation. Consequently, the :ref:`ParticleMesh`
representation consists of *explicit* :ref:`Geometry` and *implicit*
:ref:`Topology`. As with the ``mint::CurvilinearMesh`` and
``mint::UnstructuredMesh``, the explicit :ref:`Geometry` of the
:ref:`ParticleMesh` is represented by employing the ``mint::MeshCoordinates``
as an internal class member.

The following code snippet provides a simple examples illustrating how to
construct and operate on a :ref:`ParticleMesh`.

.. literalinclude:: ../../../examples/mint_particle_mesh.cpp
   :language: C++
   :linenos:

.. #############################################################################
..  MESH STORAGE
.. #############################################################################

.. _MeshStorageManagement:

Mesh Storage Management
------------------------

Mint provides a flexible :ref:`MeshStorageManagement` system that can
optionally interoperate with `Sidre`_ as the underlying, in-memory,
hierarchichal datastore. This enables Mint to natively conform to
`Conduit`_'s `Blueprint`_ protocol for representing a computational mesh in
memory and thereby, facilitate with the integration across different physics
packages.

Mint's :ref:`MeshStorageManagement` substrate supports three storage options.
The applicable operations and ownership state of each storage option are
summarized in the table below, followed by a brief description of each option.

.. |check| unicode:: U+2713

.. raw:: html

      <center>

+------------------------+---------+------------+--------------+
|                        | Modify  | Reallocate | Ownership    |
|                        |         |            |              |
+========================+=========+============+==============+
| :ref:`NativeStorage`   | |check| |  |check|   | Mint         |
+------------------------+---------+------------+--------------+
| :ref:`ExternalStorage` | |check| |            | Application  |
+------------------------+---------+------------+--------------+
| :ref:`SidreStorage`    | |check| |  |check|   |  `Sidre`_    |
+------------------------+---------+------------+--------------+

.. raw:: html

      </center>


.. _NativeStorage:

Native Storage
^^^^^^^^^^^^^^

A Mint object using :ref:`NativeStorage` owns all memory and associated data.
The data can be modified and the associated memory space can be reallocated
to grow and shrink as needed. However, once the Mint object goes out-of-scope,
all data is deleted and the memory is returned to the system.

See the :ref:`sections/tutorial` for more information and a set of concrete
examples on :ref:`constructingAMesh` using :ref:`NativeStorage`.

.. _ExternalStorage:

External Storage
^^^^^^^^^^^^^^^^

A Mint object using :ref:`ExternalStorage` has a pointer to a supplied
application buffer. In this case, the data can be modified, but, the application
maintains ownership of the underlying memory. Consequently, the memory space
cannot be reallocated and once the Mint object goes out-of-scope, the data is not
deleted. The data remains persistent in the application buffers until it is
deleted by the application.

See the :ref:`sections/tutorial` for more information and a set of concrete
examples on :ref:`constructAMeshFromExternal`.

.. _SidreStorage:

Sidre Storage
^^^^^^^^^^^^^

A Mint object using :ref:`SidreStorage` is associated with a `Sidre`_ Group
object which has owneship of the mesh data. In this case the data can be
modified and the associated memory can be reallocated to grow and shrink as
needed. However, when the Mint object goes out-of-scope, the data remains
persistent in `Sidre`_.

See the :ref:`sections/tutorial` for more information and a set of concrete
examples on :ref:`usingMintWithSidre`.

.. #############################################################################
..  FOOTNOTES
.. #############################################################################

.. rubric:: Footnotes

.. [#f1] A *Mesh* is also sometimes referred to as a *Grid*
.. [#f2] The cell numbering convention is also referred to as winding convention.
.. [#f3] Local adaptive refinement is also more broadly known as local h-refinement.

.. #############################################################################
..  CITATIONS
.. #############################################################################

.. include:: citations.rst
