// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_UNSTRUCTUREDMESH_HPP_
#define MINT_UNSTRUCTUREDMESH_HPP_

// Axom includes
#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

// Mint includes
#include "axom/mint/config.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/mesh/ConnectivityArray.hpp"
#include "axom/mint/mesh/Mesh.hpp"
#include "axom/mint/mesh/MeshCoordinates.hpp"
#include "axom/mint/mesh/MeshTypes.hpp"
#include "axom/mint/mesh/internal/MeshHelpers.hpp"

// Slic includes
#include "axom/slic/interface/slic.hpp"

// C/C++ includes
#include <cstring>  // for std::memcpy

namespace axom
{
namespace mint
{
enum Topology
{
  SINGLE_SHAPE,
  MIXED_SHAPE
};

template <Topology TOPO>
struct topology_traits
{ };

template <>
struct topology_traits<SINGLE_SHAPE>
{
  constexpr static ConnectivityType cell_to_nodes = NO_INDIRECTION;
  constexpr static ConnectivityType cell_to_faces = NO_INDIRECTION;
  constexpr static ConnectivityType face_to_nodes = NO_INDIRECTION;
};

template <>
struct topology_traits<MIXED_SHAPE>
{
  constexpr static ConnectivityType cell_to_nodes = TYPED_INDIRECTION;
  constexpr static ConnectivityType cell_to_faces = INDIRECTION;
  constexpr static ConnectivityType face_to_nodes = TYPED_INDIRECTION;
};

/*!
 * \class UnstructuredMesh
 *
 * \brief Provides the ability to store and operate on Unstructured meshes.
 *
 *  An <em> unstructured mesh </em> stores both node and topology information
 *  explicitly. This allows the flexibility of discretizing the solution
 *  domain using a variety of cell types not just quadrilateral (in 2D) or
 *  hexahedral (in 3D) cells. Due to this added flexibility, the use of
 *  <em> unstructured meshes </em> is more common when dealing with complex
 *  geometries. However, <em> unstructured meshes </em> require additional
 *  storage and generally incur some performance penalty to store, create and
 *  access mesh topology information respectively.
 *
 *  Mint classifies <em> unstructured meshes </em> in two basic types based on
 *  the underlying mesh topology:
 *
 *  * <b> Single Shape Topology </b>
 *
 *   In this case, the <em> unstructured mesh </em> consists of a single cell
 *   type, e.g., a quad or triangle mesh in 2D, or, a hex or tet mesh in 3D.
 *   In this case  the underlying implementation is optimized for the
 *   specified cell type (specified in the constructor).
 *
 *  * <b> Mixed Shape Topology </b>
 *
 *    When <em> mixed cell topology </em> is specified, the <em> unstructured
 *    mesh </em> can be composed of any of the supported cell types, e.g.,
 *    a mesh consisting of both quads and triangles. This mode incurs
 *    additional overhead for storage and access to mesh topology information,
 *    since it requires indirection.
 *
 *  The list of supported cell types for an <em> unstructured mesh </em> is
 *  available in CellTypes.hpp
 *
 * An UnstructuredMesh object may be constructed using (a) native storage,
 * (b) external storage, or, (c) with Sidre storage when Mint is compiled
 *  with Sidre support:
 *
 *  * <b> Native Storage </b> <br />
 *
 *    When using native storage, the Unstructured object owns all memory
 *    associated with the particle data. The storage can grow dynamically as
 *    needed by the application, i.e., adding more nodes/cells. Once the
 *    mesh object goes out-of-scope, all memory associated with it is
 *    returned to the system.
 *
 *  * <b> External Storage </b> <br />
 *
 *    An UnstructuredMesh may also be constructed from external, user-supplied
 *    buffers. In this case, all memory associated with the mesh data
 *    is owned by the caller. Consequently, the mesh cannot grow dynamically.
 *    All calls to resize(), append(), etc. will fail with an error.
 *
 *  * <b> Sidre </b> <br />
 *
 *    An UnstructuredMesh may also be constructed using Sidre as the back-end
 *    data-store manager. The mesh data is laid out in Sidre according to the
 *    <a href="http://llnl-conduit.readthedocs.io/en/latest/"> computational
 *    mesh blueprint </a> conventions. In this case, all operations are
 *    supported, including dynamically resizing the mesh and growing the
 *    associated storage as needed. However, Sidre owns all the memory. Once the
 *    mesh object goes out-of-scope, the data remains persistent in Sidre.
 *
 * \see mint::Mesh
 * \see mint::CellTypes
 */
template <Topology TOPO>
class UnstructuredMesh : public Mesh
{
  AXOM_STATIC_ASSERT(TOPO == SINGLE_SHAPE || TOPO == MIXED_SHAPE);

public:
  /*! \brief The types for face-cell and cell-face connectivity.
   *
   * Usually, each face will connect two cells.  The exceptions are
   * - edge faces, which will connect to one cell.  The other will be -1.
   * - some faces in a malformed mesh, which may join more than two cells.
   *   The concensus is that storing such a non-manifold mesh is not useful
   *   so we only store the first two incident cells and return an error.
   *
   * Face types are stored in the face-node connectivity.
   */
  using CellToNodeConnectivity =
    ConnectivityArray<topology_traits<TOPO>::cell_to_nodes>;
  using CellToFaceConnectivity =
    ConnectivityArray<topology_traits<TOPO>::cell_to_faces>;
  using FaceToCellConnectivity = ConnectivityArray<NO_INDIRECTION>;
  using FaceToNodeConnectivity =
    ConnectivityArray<topology_traits<TOPO>::face_to_nodes>;

  /*!
   * \brief Default constructor. Disabled.
   */
  UnstructuredMesh() = delete;

  /// \name Native Storage Constructors
  /// @{

  /*!
   * \brief Constructs an Unstructured single topology mesh.
   *
   * \param [in] ndims the number of dimensions.
   * \param [in] cell_type the cell type of the mesh.
   * \param [in] node_capacity the number of nodes to allocate space for.
   * \param [in] cell_capacity the number of cells to allocate space for.
   *
   * \note This constructor is only active when TOPO == SINGLE_SHAPE.
   *
   * \post getCellType() == cell_type
   * \post getNumberOfNodes() == 0
   * \post getNumberOfCells() == 0
   */
  UnstructuredMesh(int ndims,
                   CellType cell_type,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT)
    : Mesh(ndims, UNSTRUCTURED_MESH)
    , m_coordinates(new MeshCoordinates(ndims, 0, node_capacity))
    , m_cell_to_node(new CellToNodeConnectivity(cell_type, cell_capacity))
    , m_cell_to_face(initializeCellToFace(cell_type))
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode(cell_type))
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == SINGLE_SHAPE,
      "This constructor is only active for single topology meshes.");

    SLIC_ERROR_IF(
      cell_type == PRISM || cell_type == PYRAMID,
      "Single shape unstructured meshes do not support prisms or pyramids");

    initialize();
  }

  /*!
   * \brief Constructs an Unstructured mixed topology mesh.
   *
   * \param [in] ndims the number of dimensions.
   * \param [in] node_capacity the number of nodes to allocate space for.
   * \param [in] cell_capacity the number of cells to allocate space for.
   * \param [in] connectivity_capacity the number of vertices to allocate space
   *  for in the cell connectivity array.
   *
   * \note This constructor is only active when TOPO == MIXED_SHAPE.
   *
   * \post getCellType() == UNDEFINED_CELL
   * \post getNumberOfNodes() == 0
   * \post getNumberOfCells() == 0
   */
  UnstructuredMesh(int ndims,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT,
                   IndexType connectivity_capacity = USE_DEFAULT)
    : Mesh(ndims, UNSTRUCTURED_MESH)
    , m_coordinates(new MeshCoordinates(ndims, 0, node_capacity))
    , m_cell_to_node(
        new CellToNodeConnectivity(cell_capacity, connectivity_capacity))
    , m_cell_to_face(initializeCellToFace())
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode())
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == MIXED_SHAPE,
      "This constructor is only active for mixed topology meshes.");

    m_has_mixed_topology = true;
    initialize();
  }

  /// @}

  /// \name External Storage Constructors
  /// @{

  /*!
   * \brief Constructs an Unstructured single topology mesh using the provided
   *  external buffers.
   *
   * \param [in] cell_type the cell type of the mesh.
   * \param [in] n_cells the number of cells in the mesh.
   * \param [in] cell_capacity max number of cells the mesh is able to hold.
   * \param [in] connectivity the cell connectivity array.
   * \param [in] n_nodes the number of nodes in the mesh.
   * \param [in] node_capacity max number of nodes the mesh is able to hold
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates (required only for 2D and 3D)
   * \param [in] z pointer to the z-coordinates (required only for 3D)
   *
   * \note The length of the connectivity array must be at least
   *  cell_capacity * getCellInfo( cell_type ).num_nodes.
   * \note The provided coordinate arrays have to be of length at least
   *  node_capacity.
   *
   * \note This constructor is only supported when TOPO == SINGLE_SHAPE.
   *
   * \post getCellType() == cell_type
   * \post getNumberOfCells() == n_cells
   * \post getCellCapacity == cell_capacity
   * \post getNumberOfNodes() == n_nodes
   * \post getNodeCapacity() == node_capacity
   * \post isExternal() == true
   */
  UnstructuredMesh(CellType cell_type,
                   IndexType n_cells,
                   IndexType cell_capacity,
                   IndexType* connectivity,
                   IndexType n_nodes,
                   IndexType node_capacity,
                   double* x,
                   double* y = nullptr,
                   double* z = nullptr)
    : Mesh(internal::dim(x, y, z), UNSTRUCTURED_MESH)
    , m_coordinates(new MeshCoordinates(n_nodes, node_capacity, x, y, z))
    , m_cell_to_node(
        new CellToNodeConnectivity(cell_type, n_cells, connectivity, cell_capacity))
    , m_cell_to_face(initializeCellToFace(cell_type))
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode(cell_type))
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == SINGLE_SHAPE,
      "This constructor is only active for single topology meshes.");

    SLIC_ERROR_IF(
      cell_type == PRISM || cell_type == PYRAMID,
      "Single shape unstructured meshes do not support prisms or pyramids");

    SLIC_ASSERT(x != nullptr);
    SLIC_ASSERT(m_ndims < 2 || y != nullptr);
    SLIC_ASSERT(m_ndims < 3 || z != nullptr);

    initialize();
  }

  /*!
   * \brief Constructs an Unstructured single topology mesh using the provided
   *  external buffers.
   *
   * \param [in] cell_type the cell type of the mesh
   * \param [in] n_cells the number of cells in the mesh
   * \param [in] connectivity the cell connectivity array.
   * \param [in] n_nodes the number of nodes in the mesh
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates (required only for 2D or 3D)
   * \param [in] z pointer to the z-coordinates (required only for 3D)
   *
   * \note This constructor is only supported when TOPO == SINGLE_SHAPE.
   *
   * \post getCellType() == cell_type
   * \post getNumberOfCells() == n_cells
   * \post getCellCapacity == n_cells
   * \post getNumberOfNodes() == n_nodes
   * \post getNodeCapacity() == n_nodes
   * \post isExternal() == true
   */
  UnstructuredMesh(CellType cell_type,
                   IndexType n_cells,
                   IndexType* connectivity,
                   IndexType n_nodes,
                   double* x,
                   double* y = nullptr,
                   double* z = nullptr)
    : UnstructuredMesh(cell_type,
                       n_cells,
                       n_cells,
                       connectivity,
                       n_nodes,
                       n_nodes,
                       x,
                       y,
                       z)
  { }

  /*!
   * \brief Constructs an Unstructured mixed topology mesh using the provided
   *  external buffers.
   *
   * \param [in] n_cells the number of cells in the mesh.
   * \param [in] cell_capacity max number of cells the mesh is able to hold
   * \param [in] connectivity_capacity max capacity of the connectivity array
   * \param [in] connectivity the connectivity array.
   * \param [in] offsets array of cell offsets (of length cell_capacity + 1)
   * \param [in] types array of cell types (of length cell_capacity)
   * \param [in] n_nodes the number of nodes in the mesh.
   * \param [in] node_capacity max number of nodes the mesh can hold
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates (required only for 2D and 3D)
   * \param [in] z pointer to the z-coordinates (required only 3D)
   *
   * \note The supplied connectivity array must have a length that is at least
   *  equal to the specified connectivity_capacity.
   * \note the provided coordinate arrays are to be of length at least
   *  node_capacity.
   *
   * \note This constructor is only supported when TOPO == MIXED_SHAPE.
   *
   * \post getCellType() == UNDEFINED_CELL
   * \post getNumberOfCells() == n_cells
   * \post getCellCapacity == cell_capacity
   * \post getNumberOfNodes() == n_nodes
   * \post getNodeCapacity() == node_capacity
   * \post isExternal() == true
   */
  UnstructuredMesh(IndexType n_cells,
                   IndexType cell_capacity,
                   IndexType connectivity_capacity,
                   IndexType* connectivity,
                   IndexType* offsets,
                   CellType* types,
                   IndexType n_nodes,
                   IndexType node_capacity,
                   double* x,
                   double* y = nullptr,
                   double* z = nullptr)
    : Mesh(internal::dim(x, y, z), UNSTRUCTURED_MESH)
    , m_coordinates(new MeshCoordinates(n_nodes, node_capacity, x, y, z))
    , m_cell_to_node(new CellToNodeConnectivity(n_cells,
                                                connectivity,
                                                offsets,
                                                types,
                                                cell_capacity,
                                                connectivity_capacity))
    , m_cell_to_face(initializeCellToFace())
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode())
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == MIXED_SHAPE,
      "This constructor is only active for mixed topology meshes.");

    SLIC_ASSERT(x != nullptr);
    SLIC_ASSERT(m_ndims < 2 || y != nullptr);
    SLIC_ASSERT(m_ndims < 3 || z != nullptr);

    m_has_mixed_topology = true;
    initialize();
  }

  /*!
   * \brief Constructs an Unstructured mixed topology mesh using the provided
   *  external mesh buffers.
   *
   * \param [in] n_cells the number of cells in the mesh.
   * \param [in] connectivity_size the size of the connectivity array
   * \param [in] connectivity the connectivity array
   * \param [in] offsets array of cell offsets (of length n_cells+1)
   * \param [in] types array of cell types (of length n_cells)
   * \param [in] n_nodes the number of nodes in the mesh
   * \param [in] x pointer to the x-coordinates
   * \param [in] y pointer to the y-coordinates (required only for 2D and 3D)
   * \param [in] z pointer to the z-coordinates (required only for 3D)
   *
   * \note The supplied connectivity array must have a length that is at least
   *  equal to the specified connectivity_capacity.
   * \note the provided coordinate arrays are to be of length at least
   *  node_capacity.
   *
   * \post getCellType() == UNDEFINED_CELL
   * \post getNumberOfCells() == n_cells
   * \post getCellCapacity == n_cells
   * \post getNumberOfNodes() == n_nodes
   * \post getNodeCapacity() == n_cells
   * \post isExternal() == true
   */
  UnstructuredMesh(IndexType n_cells,
                   IndexType connectivity_size,
                   IndexType* connectivity,
                   IndexType* offsets,
                   CellType* types,
                   IndexType n_nodes,
                   double* x,
                   double* y = nullptr,
                   double* z = nullptr)
    : UnstructuredMesh(n_cells,
                       n_cells,
                       connectivity_size,
                       connectivity,
                       offsets,
                       types,
                       n_nodes,
                       n_nodes,
                       x,
                       y,
                       z)
  { }

  /// @}

  /// \name Sidre Storage Constructors
  /// @{

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Creates an UnstructuredMesh instance from a given Sidre group that
   *  holds mesh data for an unstructured mesh according to the conventions
   *  described in the computational mesh blueprint.
   *
   * \param [in] group the sidre::Group to use.
   * \param [in] topo optional argument specifying the name of the topology
   *  associated with this Mesh instance.
   *
   * \note If a topology name is not provided, the implementation will construct
   *  a mesh based on the 1st topology group under the parent "topologies"
   *  group.
   *
   * \pre group != nullptr.
   * \pre blueprint::isValidRootGroup( group ) == true
   * \post isInSidre() == true
   */
  UnstructuredMesh(sidre::Group* group, const std::string& topo = "")
    : Mesh(group, topo)
    , m_coordinates(new MeshCoordinates(getCoordsetGroup()))
    , m_cell_to_node(new CellToNodeConnectivity(getTopologyGroup()))
    , m_cell_to_face(nullptr)
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(nullptr)
  {
    SLIC_ERROR_IF(
      m_type != UNSTRUCTURED_MESH,
      "Supplied sidre::Group does not correspond to a UnstructuredMesh.");

    if(TOPO == MIXED_SHAPE)
    {
      m_has_mixed_topology = true;
      m_cell_to_face = initializeCellToFace();
      m_face_to_node = initializeFaceToNode();
    }
    else
    {
      SLIC_ERROR_IF(
        getCellType() == PRISM || getCellType() == PYRAMID,
        "Single shape unstructured meshes do not support prisms or pyramids");
      m_cell_to_face = initializeCellToFace(getCellType());
      m_face_to_node = initializeFaceToNode(getCellType());
    }

    initialize();
  }

  /*!
   * \brief Creates an UnstructuredMesh instance on an empty Sidre group.
   *
   * \param [in] ndims the number of dimensions.
   * \param [in] cell_type the cell type of the mesh.
   * \param [in] group the sidre::Group to use.
   * \param [in] topo the name of the associated topology group.
   * \param [in] coordset the name of the associated coordset group.
   * \param [in] node_capacity the number of nodes to allocate space for.
   * \param [in] cell_capacity the number of cells to allocate space for.
   * \param [in] connectivity_capacity the number of vertices to allocate space
   *  for in the cell connectivity array.
   *
   * \note If a topology and coordset name are not provided a default name is
   *  used by the implementation.
   * \note The first two constructors are only active when
   *  TOPO == SINGLE_SHAPE and the last two are active only when
   *  TOPO == MIXED_SHAPE.
   *
   * \pre group != nullptr.
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   * \post blueprint::isValidRootGroup( group )
   * \post getNumberOfNodes() == 0
   * \post getNumberOfCells() == 0
   * \post isInSidre() == true
   */
  /// @{

  UnstructuredMesh(int ndims,
                   CellType cell_type,
                   sidre::Group* group,
                   const std::string& topo,
                   const std::string& coordset,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT)
    : Mesh(ndims, UNSTRUCTURED_MESH, group, topo, coordset)
    , m_coordinates(
        new MeshCoordinates(getCoordsetGroup(), ndims, 0, node_capacity))
    , m_cell_to_node(new CellToNodeConnectivity(cell_type,
                                                getTopologyGroup(),
                                                getCoordsetName(),
                                                cell_capacity))
    , m_cell_to_face(initializeCellToFace(cell_type))
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode(cell_type))
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == SINGLE_SHAPE,
      "This constructor is only active for single topology meshes.");

    SLIC_ERROR_IF(
      cell_type == PRISM || cell_type == PYRAMID,
      "Single shape unstructured meshes do not support prisms or pyramids");

    initialize();
  }

  UnstructuredMesh(int ndims,
                   CellType cell_type,
                   sidre::Group* group,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT)
    : UnstructuredMesh(ndims, cell_type, group, "", "", node_capacity, cell_capacity)
  { }

  UnstructuredMesh(int ndims,
                   sidre::Group* group,
                   const std::string& topo,
                   const std::string& coordset,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT,
                   IndexType connectivity_capacity = USE_DEFAULT)
    : Mesh(ndims, UNSTRUCTURED_MESH, group, topo, coordset)
    , m_coordinates(
        new MeshCoordinates(getCoordsetGroup(), ndims, 0, node_capacity))
    , m_cell_to_node(new CellToNodeConnectivity(getTopologyGroup(),
                                                getCoordsetName(),
                                                cell_capacity,
                                                connectivity_capacity))
    , m_cell_to_face(initializeCellToFace())
    , m_face_to_cell(new FaceToCellConnectivity(2, 0))
    , m_face_to_node(initializeFaceToNode())
  {
    AXOM_STATIC_ASSERT_MSG(
      TOPO == MIXED_SHAPE,
      "This constructor is only active for mixed topology meshes.");

    m_has_mixed_topology = true;
    initialize();
  }

  UnstructuredMesh(int ndims,
                   sidre::Group* group,
                   IndexType node_capacity = USE_DEFAULT,
                   IndexType cell_capacity = USE_DEFAULT,
                   IndexType connectivity_capacity = USE_DEFAULT)
    : UnstructuredMesh(ndims,
                       group,
                       "",
                       "",
                       node_capacity,
                       cell_capacity,
                       connectivity_capacity)
  { }

    /// @}

#endif /* AXOM_MINT_USE_SIDRE */

  /// @}

  /// \name Virtual methods
  /// @{

  /*!
   * \brief Destructor, deletes the MeshCoordinates and ConnectivityArray.
   */
  virtual ~UnstructuredMesh()
  {
    delete m_coordinates;
    m_coordinates = nullptr;

    delete m_cell_to_node;
    m_cell_to_node = nullptr;

    delete m_cell_to_face;
    m_cell_to_face = nullptr;

    delete m_face_to_cell;
    m_face_to_cell = nullptr;

    delete m_face_to_node;
    m_face_to_node = nullptr;
  }

  /// \name Cells
  /// @{

  /*!
   * \brief Return the number of cells in the mesh.
   */
  virtual IndexType getNumberOfCells() const final override
  {
    return m_cell_to_node->getNumberOfIDs();
  }

  /*!
   * \brief Return the capacity for cells.
   */
  virtual IndexType getCellCapacity() const final override
  {
    return m_cell_to_node->getIDCapacity();
  }

  /*!
   * \brief Return the type of the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored if TOPO == SINGLE_SHAPE. If TOPO == MIXED_SHAPE and no
   *  cellID is provided the returned type is UNDEFINED_CELL.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual CellType getCellType(IndexType cellID = -1) const final override
  {
    return m_cell_to_node->getIDType(cellID);
  }

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored if TOPO == SINGLE_SHAPE.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getNumberOfCellNodes(IndexType cellID = 0) const final override
  {
    return m_cell_to_node->getNumberOfValuesForID(cellID);
  }

  /*!
   * \brief Copy the connectivity of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] nodes the buffer into which the connectivity is copied.
   *
   * \return The number of nodes for the given cell.
   *
   * \pre nodes != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellNodeIDs(IndexType cellID,
                                   IndexType* nodes) const final override
  {
    SLIC_ASSERT(nodes != nullptr);
    const IndexType n_nodes = getNumberOfCellNodes(cellID);
    std::memcpy(nodes, getCellNodeIDs(cellID), n_nodes * sizeof(IndexType));
    return n_nodes;
  }

  /*!
   * \brief Return the number of faces associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   */
  virtual IndexType getNumberOfCellFaces(IndexType cellID = 0) const final override
  {
    return getCellInfo(getCellType(cellID)).num_faces;
  }

  /*!
   * \brief Copy the face IDs of the given cell into the provided buffer.
   * The buffer must be of length at least getNumberOfCellFaces( cellID ).
   *
   * \param [in] cellID the ID of the cell in question
   * \param [out] faces the buffer into which the face IDs are copied.
   *
   * \return The number of faces for the given cell.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   *
   * \pre faces != nullptr
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCellFaceIDs(IndexType cellID,
                                   IndexType* faces) const final override
  {
    SLIC_ASSERT(faces != nullptr);
    const IndexType n_faces = getNumberOfCellFaces(cellID);
    std::memcpy(faces, getCellFaceIDs(cellID), n_faces * sizeof(IndexType));
    return n_faces;
  }

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the number of nodes in the mesh.
   */
  virtual IndexType getNumberOfNodes() const final override
  {
    return m_coordinates->numNodes();
  }

  /*!
   * \brief Return the capacity for nodes.
   */
  virtual IndexType getNodeCapacity() const final override
  {
    return m_coordinates->capacity();
  }

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != nullptr
   */
  virtual void getNode(IndexType nodeID, double* coords) const final override
  {
    m_coordinates->getCoordinates(nodeID, coords);
  }

  /*!
   * \brief Return a pointer to the array of nodal coordinates of the
   *  given dimension.
   *
   * \param [in] dim the dimension to return.
   *
   * \pre 0 <= dim < getDimension()
   */
  /// @{

  virtual double* getCoordinateArray(int dim) final override
  {
    return m_coordinates->getCoordinateArray(dim);
  }

  virtual const double* getCoordinateArray(int dim) const final override
  {
    return m_coordinates->getCoordinateArray(dim);
  }

  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Return the number of faces in the mesh.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   */
  virtual IndexType getNumberOfFaces() const final override
  {
    return m_face_to_cell->getNumberOfIDs();
  }

  /*!
   * \brief Return the capacity for faces.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   */
  virtual IndexType getFaceCapacity() const final override
  {
    return m_face_to_cell->getIDCapacity();
  }

  /*!
   * \brief Return the type of the given face.
   *
   * \param [in] faceID the ID of the face in question.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   */
  virtual CellType getFaceType(IndexType faceID) const final override
  {
    return m_face_to_node->getIDType(faceID);
  }

  /*!
   * \brief Return the number of nodes associated with the given face.
   *
   * \param [in] faceID the ID of the face in question.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   */
  virtual IndexType getNumberOfFaceNodes(IndexType faceID = 0) const final override
  {
    return m_face_to_node->getNumberOfValuesForID(faceID);
  }

  /*!
   * \brief Copy the IDs of the nodes that compose the given face into the
   *  provided buffer.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] nodes the buffer into which the node IDs are copied, must
   *  be of length at least getNumberOfFaceNodes().
   *
   * \return The number of nodes for the given face.
   *
   * A face with ID faceID will have a normal pointing outward from the
   * first cell it is adjacent to, as returned by
   * getFaceCellIDs(faceID, cellID1, cellID2).
   *
   * \pre nodes != nullptr
   * \pre 0 <= faceID < getNumberOfCells()
   */
  virtual IndexType getFaceNodeIDs(IndexType faceID,
                                   IndexType* nodes) const final override
  {
    SLIC_ASSERT(nodes != nullptr);
    const IndexType n_nodes = getNumberOfFaceNodes(faceID);
    std::memcpy(nodes, getFaceNodeIDs(faceID), n_nodes * sizeof(IndexType));
    return n_nodes;
  }

  /*!
   * \brief Copy the IDs of the cells adjacent to the given face into the
   *  provided indices.
   *
   * \param [in] faceID the ID of the face in question.
   * \param [out] cellIDOne the ID of the first cell.
   * \param [out] cellIDTwo the ID of the second cell.
   *
   * \note A face can be associated with one or two cells, depending on whether
   *  it is an external boundary face or interior face. By convention, if a face
   *  is an external boundary face, then only cellIDOne exists and cellIDTwo
   *  will be set to -1.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   *
   * \pre 0 <= faceID < getNumberOfFaces()
   */
  virtual void getFaceCellIDs(IndexType faceID,
                              IndexType& cellIDOne,
                              IndexType& cellIDTwo) const final override
  {
    IndexType* faces = (*m_face_to_cell)[faceID];
    cellIDOne = faces[0];
    cellIDTwo = faces[1];
  }

  /// @}

  /// \name Edges
  /// @{

  /*!
   * \brief Return the number of edges in the mesh.
   */
  virtual IndexType getNumberOfEdges() const final override
  {
    SLIC_ERROR("NOT IMPLEMENTED!!!");
    return 0;
  }

  /*!
   * \brief Return the capacity for edges.
   */
  virtual IndexType getEdgeCapacity() const final override
  {
    SLIC_ERROR("NOT IMPLEMENTED!!!");
    return 0;
  }

  /// @}

  /*!
   * \brief Return true iff both the connectivity and coordinates are stored in
   *  external arrays.
   */
  virtual bool isExternal() const final override
  {
    bool connec_external = m_cell_to_node->isExternal();
    bool coords_external = m_coordinates->isExternal();

    if(connec_external != coords_external)
    {
      SLIC_WARNING("External state not consistent.");
      return false;
    }

    return connec_external;
  }

  /// @}

  /// \name Attribute get/set Methods
  /// @{

  /// \name Cells
  /// @{

  /*!
   * \brief Return the cell resize ratio.
   */
  double getCellResizeRatio() const { return m_cell_to_node->getResizeRatio(); }

  /*!
   * \brief Set the cell resize ratio.
   *
   * \param [in] ratio the new cell resize ratio.
   *
   * \post getCellResizeRatio() == ratio
   */
  void setCellResizeRatio(double ratio)
  {
    m_cell_to_node->setResizeRatio(ratio);
    m_mesh_fields[CELL_CENTERED]->setResizeRatio(ratio);
  }

  /*!
   * \brief Return the size of the connectivity array.
   */
  IndexType getCellNodesSize() const
  {
    return m_cell_to_node->getNumberOfValues();
  }

  /*!
   * \brief Return the capacity of the connectivity array.
   */
  IndexType getCellNodesCapacity() const
  {
    return m_cell_to_node->getValueCapacity();
  }

  /*!
   * \brief Resizes the cell connectivity array and cell-centered fields of this
   *  mesh instance to hold the specified number of cells.
   *
   * \param [in] cell_size the number of cells to resize to.
   *
   * \post getNumberOfCells() == cell_size
   */
  void resizeCells(IndexType cell_size)
  {
    IndexType connectivity_size =
      (hasMixedCellTypes()) ? USE_DEFAULT : getNumberOfCellNodes() * cell_size;
    m_cell_to_node->resize(cell_size, connectivity_size);
    m_mesh_fields[CELL_CENTERED]->resize(cell_size);
  }

  /*!
   * \brief Reserve space for the given number of cells.
   *
   * \param [in] cell_capacity the number of cells to reserve space for.
   * \param [in] connectivity_capacity the ammount of space to reserve in the
   *  connectivity array. Ignored if TOPO == SINGLE_SHAPE.
   *
   * \post getCellCapacity() >= cell_capacity
   */
  void reserveCells(IndexType cell_capacity,
                    IndexType connectivity_capacity = USE_DEFAULT)
  {
    m_cell_to_node->reserve(cell_capacity, connectivity_capacity);
    m_mesh_fields[CELL_CENTERED]->reserve(cell_capacity);
  }

  /*!
   * \brief Shrink the cell capacity to be equal to the number of cells.
   *
   * \post getCellCapacity() == getNumberOfCells()
   * \post getCellNodesCapacity() == getCellNodesSize()
   */
  void shrinkCells()
  {
    m_cell_to_node->shrink();
    m_mesh_fields[CELL_CENTERED]->shrink();
  }

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the node resize ratio.
   */
  double getNodeResizeRatio() const { return m_coordinates->getResizeRatio(); }

  /*!
   * \brief Set the node resize ratio.
   *
   * \param [in] ratio the new node resize ratio.
   *
   * \post getNodeResizeRatio() == ratio
   */
  void setNodeResizeRatio(double ratio)
  {
    m_coordinates->setResizeRatio(ratio);
    m_mesh_fields[NODE_CENTERED]->setResizeRatio(ratio);
  }

  /*!
   * \brief Resizes the nodal coordinates and fields of this mesh instance to
   *  the specified number of nodes.
   *
   * \param [in] nodes_size the number of nodes to resize to.
   *
   * \post getNumberOfNodes() == nodes_size
   */
  void resizeNodes(IndexType nodes_size)
  {
    m_coordinates->resize(nodes_size);
    m_mesh_fields[NODE_CENTERED]->resize(nodes_size);
  }

  /*!
   * \brief Reserve space for the given number of nodes.
   *
   * \param [in] node_capacity the number of nodes to reserve space for.
   *
   * \post getNodeCapacity() >= node_capacity
   */
  void reserveNodes(IndexType node_capacity)
  {
    m_coordinates->reserve(node_capacity);
    m_mesh_fields[NODE_CENTERED]->reserve(node_capacity);
  }

  /*!
   * \brief Shrink the node capacity to be equal to the number of nodes.
   *
   * \post getNodeCapacity() == getNumberOfNodes()
   */
  void shrinkNodes()
  {
    m_coordinates->shrink();
    m_mesh_fields[NODE_CENTERED]->shrink();
  }

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Return the face resize ratio.
   */
  double getFaceResizeRatio() const
  {
    const double ratio = m_face_to_node->getResizeRatio();
    SLIC_WARNING_IF(m_face_to_cell->getResizeRatio() != ratio,
                    "Resize ratios are inconsistent");
    return ratio;
  }

  /*!
   * \brief Return the size of the connectivity array.
   */
  IndexType getFaceNodesSize() const
  {
    return m_face_to_node->getNumberOfValues();
  }

  /*!
   * \brief Return the capacity of the connectivity array.
   */
  IndexType getFaceNodesCapacity() const
  {
    return m_face_to_node->getValueCapacity();
  }

  /// @}

  /// \name Edges
  /// @{

  /*!
   * \brief Return the edge resize ratio.
   */
  double getEdgeResizeRatio() const
  {
    SLIC_ERROR("NOT IMPLEMENTED!!!");
    return 0.0;
  }

  /// @}

  /*!
   * \brief Resizes this mesh instance to the specified number of nodes & cells.
   *
   * \param [in] node_size the desired number of nodes
   * \param [in] cell_size the desired number of cells
   *
   * \note This method will also resize the node-centered and cell-centered
   *  fields accordingly.
   *
   * \post getNumberOfNodes() == nodes_size
   * \post getNumberOfCells() == cell_size
   *
   * \see resizeNodes()
   * \see resizeCells()
   */
  void resize(IndexType node_size, IndexType cell_size)
  {
    resizeNodes(node_size);
    resizeCells(cell_size);
  }

  /*!
   * \brief Reserve space for the given number of nodes and cells.
   *
   * \param [in] node_capacity the number of nodes to reserve space for.
   * \param [in] cell_capacity the number of cells to reserve space for.
   * \param [in] connectivity_capacity the ammount of space to reserve in the
   *  connectivity array. Ignored if TOPO == SINGLE_SHAPE.
   *
   * \post getNodeCapacity() >= node_capacity
   * \post getCellCapacity() >= cell_capacity
   */
  void reserve(IndexType node_capacity,
               IndexType cell_capacity,
               IndexType connectivity_capacity = USE_DEFAULT)
  {
    reserveNodes(node_capacity);
    reserveCells(cell_capacity, connectivity_capacity);
  }

  /*!
   * \brief Shrink the node capacity to be equal to the number of nodes and the
   *  cell capacity to be equal to the number of cells.
   *
   * \post getNodeCapacity() == getNumberOfNodes()
   * \post getCellCapacity() == getNumberOfCells()
   * \post getCellNodesCapacity() == getCellNodesSize()
   */
  void shrink()
  {
    shrinkNodes();
    shrinkCells();
  }

  /*!
   * \brief Return true iff the mesh holds no nodes and no cells.
   */
  bool empty() const
  {
    return m_coordinates->empty() && m_cell_to_node->empty();
  }

  /*!
   * \brief Return true iff both the connectivity and coordinates are stored in
   *  sidre.
   */
  bool isInSidre() const
  {
    bool connec_sidre = m_cell_to_node->isInSidre();
    bool coords_sidre = m_coordinates->isInSidre();

    if(connec_sidre != coords_sidre)
    {
      SLIC_WARNING("Sidre state not consistent.");
      return false;
    }

    return connec_sidre;
  }

  /// @}

  /// \name Data Access Methods
  /// @{

  /// \name Cells
  /// @{

  /*!
   * \brief Return a pointer to the connectivity of the given cell. The
   *  buffer is guarenteed to be of length at least
   *  getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  /// @{

  IndexType* getCellNodeIDs(IndexType cellID)
  {
    return (*m_cell_to_node)[cellID];
  }

  const IndexType* getCellNodeIDs(IndexType cellID) const
  {
    return (*m_cell_to_node)[cellID];
  }

  /*!
   * \brief Return a pointer to the faces of the given cell. The
   *  buffer is guarenteed to be of length at least
   *  getNumberOfCellFaces( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  /// @{

  IndexType* getCellFaceIDs(IndexType cellID)
  {
    return (*m_cell_to_face)[cellID];
  }

  const IndexType* getCellFaceIDs(IndexType cellID) const
  {
    return (*m_cell_to_face)[cellID];
  }

  /// @}

  /*!
   * \brief Return a pointer to the cell nodes array, of length
   *  getCellNodesSize().
   */
  /// @{

  IndexType* getCellNodesArray() { return m_cell_to_node->getValuePtr(); }

  const IndexType* getCellNodesArray() const
  {
    return m_cell_to_node->getValuePtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the cell nodes offset array, of length
   *  getNumberOfCells() + 1. Returns nullptr if
   *  TOPO == SINGLE_SHAPE.
   */
  /// @{

  IndexType* getCellNodesOffsetsArray()
  {
    return m_cell_to_node->getOffsetPtr();
  }

  const IndexType* getCellNodesOffsetsArray() const
  {
    return m_cell_to_node->getOffsetPtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the cell types array, of length
   *  getNumberOfCells(). Returns nullptr if
   *  TOPO == SINGLE_SHAPE.
   */
  /// @{

  CellType* getCellTypesArray() { return m_cell_to_node->getTypePtr(); }

  const CellType* getCellTypesArray() const
  {
    return m_cell_to_node->getTypePtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the cell faces array, of length
   *  getCellNodesSize().
   */
  /// @{

  IndexType* getCellFacesArray() { return m_cell_to_face->getValuePtr(); }

  const IndexType* getCellFacesArray() const
  {
    return m_cell_to_face->getValuePtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the cell faces offset array, of length
   *  getNumberOfCells() + 1.
   */
  /// @{

  IndexType* getCellFacesOffsetsArray()
  {
    return m_cell_to_face->getOffsetPtr();
  }

  const IndexType* getCellFacesOffsetsArray() const
  {
    return m_cell_to_face->getOffsetPtr();
  }

  /// @}

  /*!
   * \brief Append a cell to the mesh.
   *
   * \param [in] connec the connectivity of the new cell.
   * \param [in] type the type of the new cell, ignored if
   *  TOPO == SINGLE_SHAPE.
   *
   * \pre connec != nullptr
   */
  void appendCell(const IndexType* connec, CellType type = UNDEFINED_CELL)
  {
    IndexType n_values =
      (type == UNDEFINED_CELL) ? 0 : getCellInfo(type).num_nodes;
    m_cell_to_node->append(connec, n_values, type);
    m_mesh_fields[CELL_CENTERED]->resize(getNumberOfCells());
  }

  /*!
   * \brief Append multiple cells to the mesh.
   *
   * \param [in] connec the connectivity of the new cells.
   * \param [in] n_cells the number of cells to append.
   * \param [in] offsets the offsets array of the cells to append, ignored
   *  if TOPO == SINGLE_SHAPE.
   * \param [in] types the types array of the new cells, ignored if
   *  TOPO == SINGLE_SHAPE.
   *
   * \pre connec != nullptr
   * \pre n_cells >= 0
   */
  void appendCells(const IndexType* connec,
                   IndexType n_cells,
                   const IndexType* offsets = nullptr,
                   const CellType* types = nullptr)
  {
    m_cell_to_node->appendM(connec, n_cells, offsets, types);
    m_mesh_fields[CELL_CENTERED]->resize(getNumberOfCells());
  }

  /*!
   * \brief Insert a cell in to the mesh at the given position.
   *
   * \param [in] connec the connectivity of the new cell.
   * \param [in] ID the position to insert at.
   * \param [in] n_values the number of values in the connectivity, ignored
   *  if TOPO == SINGLE_SHAPE.
   * \param [in] type the type of the new cells, ignored if
   *  TOPO == SINGLE_SHAPE.
   *
   * \pre connec != nullptr
   * \pre 0 <= ID <= getNumberOfCells()
   */
  void insertCell(const IndexType* connec,
                  IndexType ID,
                  CellType type = UNDEFINED_CELL)
  {
    IndexType n_values =
      (type == UNDEFINED_CELL) ? 0 : getCellInfo(type).num_nodes;
    m_cell_to_node->insert(connec, ID, n_values, type);
    m_mesh_fields[CELL_CENTERED]->emplace(ID, 1);
  }

  /*!
   * \brief Insert multiple cells in to the mesh at the given position.
   *
   * \param [in] connec the connectivity of the new cells.
   * \param [in] start_ID the position to insert at.
   * \param [in] n_cells the number of cells to insert
   * \param [in] offsets the offsets array of the cells to append, ignored
   *  if TOPO == SINGLE_SHAPE.
   * \param [in] types the types array of the new cells, ignored if
   *  TOPO == SINGLE_SHAPE.
   *
   * \pre connec != nullptr
   * \pre 0 <= start_ID <= getNumberOfCells()
   */
  void insertCells(const IndexType* connec,
                   IndexType start_ID,
                   IndexType n_cells,
                   const IndexType* offsets = nullptr,
                   const CellType* types = nullptr)
  {
    m_cell_to_node->insertM(connec, start_ID, n_cells, offsets, types);
    m_mesh_fields[CELL_CENTERED]->emplace(start_ID, n_cells);
  }

  /// @}

  /// \name Nodes
  /// @{

  /*!
   * \brief Return the coordinate of the given dimension of the given node.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] dim the dimension to return.
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre 0 <= dim < getDimension()
   */
  double getNodeCoordinate(IndexType nodeID, int dim) const
  {
    return m_coordinates->getCoordinate(nodeID, dim);
  }

  /*!
   * \brief Appends a new node to the mesh.
   *
   * \param [in] x the first coordinate to append.
   * \param [in] y the second coordinate to append.
   * \param [in] z the third coordinate to append.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   */
  /// @{

  IndexType appendNode(double x)
  {
    IndexType n_index = m_coordinates->append(x);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
    return n_index;
  }

  IndexType appendNode(double x, double y)
  {
    IndexType n_index = m_coordinates->append(x, y);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
    return n_index;
  }

  IndexType appendNode(double x, double y, double z)
  {
    IndexType n_index = m_coordinates->append(x, y, z);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
    return n_index;
  }

  /// @}

  /*!
   * \brief Appends multiple nodes to the mesh.
   *
   * \param [in] coords pointer to the nodes to append, of length
   *  n * getDimension().
   * \param [in] n the number of nodes to append.
   *
   * \note coords is assumed to be in the array of structs format, ie
   *  coords = {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}.
   *
   * \pre coords != nullptr
   * \pre n >= 0
   */
  void appendNodes(const double* coords, IndexType n = 1)
  {
    m_coordinates->append(coords, n);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
  }

  /*!
   * \brief Appends new nodes to the mesh.
   *
   * \param [in] x array of the first coordinates to append, of length n.
   * \param [in] y array of the second coordinates to append, of length n.
   * \param [in] z array of the third coordinates to append, of length n.
   * \param [in] n the number of coordinates to append.
   *
   * \note The first method is only valid for 2D meshes while the second
   *  is only for 3D.
   * \pre x != nullptr
   * \pre y != nullptr
   * \pre z != nullptr
   * \pre n >= 0
   */
  /// @{

  void appendNodes(const double* x, const double* y, IndexType n)
  {
    m_coordinates->append(x, y, n);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
  }

  void appendNodes(const double* x, const double* y, const double* z, IndexType n)
  {
    m_coordinates->append(x, y, z, n);
    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
  }

  /// @}

  /*!
   * \brief Insert a node to the mesh.
   *
   * \param [in] nodeID the position to insert at.
   * \param [in] x the value of the first coordinate to insert.
   * \param [in] y the value of the second coordinate to insert.
   * \param [in] z the value of the third coordinate to insert.
   * \param [in] update_connectivity if true will update the connectivity so
   *  that all elements remain connected to the same coordinates as before.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   * \pre 0 <= nodeID <= getNumberOfNodes
   */
  /// @{

  void insertNode(IndexType nodeID, double x, bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, x);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, 1);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, 1);
    }
  }

  void insertNode(IndexType nodeID,
                  double x,
                  double y,
                  bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, x, y);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, 1);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, 1);
    }
  }

  void insertNode(IndexType nodeID,
                  double x,
                  double y,
                  double z,
                  bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, x, y, z);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, 1);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, 1);
    }
  }

  /// @}

  /*!
   * \brief Inserts multiple nodes to the mesh.
   *
   * \param [in] coords pointer to the nodes to insert, of length
   *  n * getDimension().
   * \param [in] n the number of nodes to append.
   * \param [in] update_connectivity if true will update the connectivity so
   *  that all elements remain connected to the same coordinates as before.
   *
   * \note coords is assumed to be in the array of structs format, ie
   *  coords = {x0, y0, z0, x1, y1, z1, ..., xn, yn, zn}.
   *
   * \pre 0 <= nodeID <= getNumberOfNodes
   * \pre coords != nullptr
   * \pre n >= 0
   */
  void insertNodes(IndexType nodeID,
                   const double* coords,
                   IndexType n = 1,
                   bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, coords, n);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, n);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, n);
    }
  }

  /*!
   * \brief Insert multiple nodes to the mesh.
   *
   * \param [in] nodeID the position to insert at.
   * \param [in] x the array of the first coordinates to insert.
   * \param [in] y the array of the second coordinates to insert.
   * \param [in] z the array of the third coordinates to insert.
   * \param [in] n the number of nodes to insert.
   * \param [in] update_connectivity if true will update the connectivity so
   *  that all elements remain connected to the same coordinates as before.
   *
   * \note The first method is only valid for 2D meshes while the second
   *  is only for 3D.
   * \pre 0 <= nodeID <= getNumberOfNodes
   * \pre x != nullptr
   * \pre y != nullptr if 2-D or 3-D
   * \pre z != nullptr if 3-D
   * \pre n >= 0
   */
  /// @{

  void insertNodes(IndexType nodeID,
                   const double* x,
                   const double* y,
                   IndexType n,
                   bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, x, y, n);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, n);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, n);
    }
  }

  void insertNodes(IndexType nodeID,
                   const double* x,
                   const double* y,
                   const double* z,
                   IndexType n,
                   bool update_connectivity = true)
  {
    m_coordinates->insert(nodeID, x, y, z, n);
    m_mesh_fields[NODE_CENTERED]->emplace(nodeID, n);
    if(update_connectivity)
    {
      cellConnectivityUpdateInsert(nodeID, n);
    }
  }

  /// @}

  /// @}

  /// \name Faces
  /// @{

  /*!
   * \brief Sets up cell-face, face-cell, and face-node connectivity.
   *
   * \param [in] force re-initialize face-related connectivity, even if it
   *             has already been done.
   */
  bool initializeFaceConnectivity(bool force = false)
  {
    if(getDimension() == 1)
    {
      return true;
    }

    if(!force && getNumberOfFaces() > 0)
    {
      return true;
    }

    if(force)
    {
      m_cell_to_face->resize(0, 0);
      m_face_to_cell->resize(0, 0);
      m_face_to_node->resize(0, 0);
    }

    IndexType facecount = 0;
    IndexType* f2cdata = nullptr;
    IndexType* c2fdata = nullptr;
    IndexType* c2ndata = nullptr;
    IndexType* c2foffsets = nullptr;
    IndexType* f2ndata = nullptr;
    IndexType* f2noffsets = nullptr;
    CellType* f2ntypes = nullptr;

    bool retval = internal::initFaces(this,
                                      facecount,
                                      f2cdata,
                                      c2fdata,
                                      c2ndata,
                                      c2foffsets,
                                      f2ndata,
                                      f2noffsets,
                                      f2ntypes);

    if(retval)
    {
      // Copy in the face connectivity data.
      const IndexType cellCount = getNumberOfCells();
      m_cell_to_face->reserve(cellCount, c2foffsets[cellCount]);
      m_cell_to_face->appendM(c2fdata, cellCount, c2foffsets);

      m_face_to_cell->reserve(facecount);
      m_face_to_cell->appendM(f2cdata, facecount);

      m_face_to_node->reserve(facecount, f2noffsets[facecount]);
      m_face_to_node->appendM(f2ndata, facecount, f2noffsets, f2ntypes);
    }

    m_mesh_fields[FACE_CENTERED]->resize(getNumberOfFaces());

    delete[] f2cdata;
    delete[] c2fdata;
    delete[] c2ndata;
    delete[] c2foffsets;
    delete[] f2ndata;
    delete[] f2noffsets;
    delete[] f2ntypes;

    return retval;
  }

  /*!
   * \brief Return a pointer to the nodes of the given face. The
   *  buffer is guaranteed to be of length at least
   *  getNumberOfFaceNodes( faceID ).
   *
   * \param [in] faceID the ID of the face in question.
   *
   * \note Codes must call initializeFaceConnectivity() before calling
   *       this method.
   *
   * \pre 0 <= faceID < getNumberOfFaces()
   */
  /// @{

  IndexType* getFaceNodeIDs(IndexType faceID)
  {
    return (*m_face_to_node)[faceID];
  }

  const IndexType* getFaceNodeIDs(IndexType faceID) const
  {
    return (*m_face_to_node)[faceID];
  }

  /// @}

  /*!
   * \brief Return a pointer to the face nodes array, of length
   *  getFaceNodesSize().
   */
  /// @{

  IndexType* getFaceNodesArray() { return m_face_to_node->getValuePtr(); }

  const IndexType* getFaceNodesArray() const
  {
    return m_face_to_node->getValuePtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the face nodes offset array, of length
   *  getNumberOfFaces() + 1.
   */
  /// @{

  IndexType* getFaceNodesOffsetsArray()
  {
    return m_face_to_node->getOffsetPtr();
  }

  const IndexType* getFaceNodesOffsetsArray() const
  {
    return m_face_to_node->getOffsetPtr();
  }

  /// @}

  /*!
   * \brief Return a pointer to the face cells array, of length
   *  2 * getNumberOfFaces().
   */
  /// @{

  IndexType* getFaceCellsArray() { return m_face_to_cell->getValuePtr(); }

  const IndexType* getFaceCellsArray() const
  {
    return m_face_to_cell->getValuePtr();
  }

  /// @}

  /// @}

  /// @}

private:
  /*! \brief Construct and fill the cell-to-face connectivity. */
  void buildCellFaceConnectivity(IndexType* c2fdata, IndexType* c2foffsets);

  /*! \brief Return a new empty CellToFaceConnectivty instance. */
  CellToFaceConnectivity* initializeCellToFace(
    CellType cell_type = UNDEFINED_CELL) const;

  /*! \brief Return a new empty FaceToNodeConnectivity instance. */
  FaceToNodeConnectivity* initializeFaceToNode(
    CellType cell_type = UNDEFINED_CELL) const;

  /*!
   * \brief Update the connectivity given an nodal insert at position pos of
   *  length n.
   *
   * \param [in] pos the position of the insert.
   * \param [in] n the length of the insert.
   */
  void cellConnectivityUpdateInsert(IndexType pos, IndexType n)
  {
    SLIC_ASSERT(0 <= pos && pos < getNumberOfNodes());

    const IndexType n_values = getCellNodesSize();
    IndexType* values = getCellNodesArray();
    SLIC_ASSERT(n_values == 0 || values != nullptr);

    for(IndexType i = 0; i < n_values; ++i)
    {
      if(values[i] >= pos)
      {
        values[i] += n;
      }
    }
  }

  /*!
   * \brief Performs common initialization.
   */
  void initialize()
  {
    m_explicit_coords = true;
    m_explicit_connectivity = true;
    m_mesh_fields[NODE_CENTERED]->setResizeRatio(getNodeResizeRatio());
    m_mesh_fields[CELL_CENTERED]->setResizeRatio(getCellResizeRatio());
    m_mesh_fields[FACE_CENTERED]->setResizeRatio(getFaceResizeRatio());

    m_mesh_fields[NODE_CENTERED]->reserve(getNodeCapacity());
    m_mesh_fields[CELL_CENTERED]->reserve(getCellCapacity());
    m_mesh_fields[FACE_CENTERED]->reserve(getFaceCapacity());

    m_mesh_fields[NODE_CENTERED]->resize(getNumberOfNodes());
    m_mesh_fields[CELL_CENTERED]->resize(getNumberOfCells());
    m_mesh_fields[FACE_CENTERED]->resize(getNumberOfFaces());
  }

  MeshCoordinates* m_coordinates;

  /*! \brief The nodes for each cell */
  CellToNodeConnectivity* m_cell_to_node;

  /*! \brief Each cell's faces */
  CellToFaceConnectivity* m_cell_to_face;

  /*! \brief Each face's cells */
  FaceToCellConnectivity* m_face_to_cell;

  /*! \brief Each face's nodes */
  FaceToNodeConnectivity* m_face_to_node;

  DISABLE_COPY_AND_ASSIGNMENT(UnstructuredMesh);
  DISABLE_MOVE_AND_ASSIGNMENT(UnstructuredMesh);
};

/*!
 * \brief Return a new CellToFaceConnectivity. This instance is for
 *  SINGLE_SHAPE meshes. If the mesh has cell type VERTEX or SEGMENT then the
 *  stride of the returned connectivity is 1 instead of 0 which is the proper
 *  number of cell faces, this avoids an error with 0 stride.
 *
 * \param [in] cell_type the cell_type of the mesh.
 */
template <>
inline UnstructuredMesh<SINGLE_SHAPE>::CellToFaceConnectivity*
UnstructuredMesh<SINGLE_SHAPE>::initializeCellToFace(CellType cell_type) const
{
  IndexType num_faces = getCellInfo(cell_type).num_faces;
  if(num_faces == 0)
  {
    num_faces = 1;
  }

  return new UnstructuredMesh<SINGLE_SHAPE>::CellToFaceConnectivity(num_faces, 0);
}

/*!
 * \brief Return a new CellToFaceConnectivity. This instance is for
 *  MIXED_SHAPE meshes.
 */
template <>
inline UnstructuredMesh<MIXED_SHAPE>::CellToFaceConnectivity*
  UnstructuredMesh<MIXED_SHAPE>::initializeCellToFace(CellType) const
{
  return new UnstructuredMesh<MIXED_SHAPE>::CellToFaceConnectivity(UNDEFINED_CELL,
                                                                   0,
                                                                   0);
}

/*!
 * \brief Return a new FaceToNodeConnectivity. This instance is for
 *  SINGLE_SHAPE meshes. If the mesh has cell type VERTEX or SEGMENT then the
 *  cell type of the returned ConnectivityArray is UNDEFINED_CELL with stride 1.
 *
 * \param [in] cell_type the cell_type of the mesh.
 */
template <>
inline UnstructuredMesh<SINGLE_SHAPE>::FaceToNodeConnectivity*
UnstructuredMesh<SINGLE_SHAPE>::initializeFaceToNode(CellType cell_type) const
{
  const CellType face_type = getCellInfo(cell_type).face_types[0];
  if(face_type == UNDEFINED_CELL)
  {
    return new UnstructuredMesh<SINGLE_SHAPE>::FaceToNodeConnectivity(1, 0);
  }

  return new UnstructuredMesh<SINGLE_SHAPE>::FaceToNodeConnectivity(face_type, 0);
}

/*!
 * \brief Return a new FaceToNodeConnectivity. This instance is for
 *  MIXED_SHAPE meshes.
 */
template <>
inline UnstructuredMesh<MIXED_SHAPE>::FaceToNodeConnectivity*
  UnstructuredMesh<MIXED_SHAPE>::initializeFaceToNode(CellType) const
{
  return new UnstructuredMesh<MIXED_SHAPE>::FaceToNodeConnectivity(0, 0);
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_UNSTRUCTUREDMESH_HPP_ */
