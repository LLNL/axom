/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MINT_MESH_HPP_
#define MINT_MESH_HPP_

#include "axom/Macros.hpp"             // for Axom macros

#include "mint/CellTypes.hpp"          // for CellType enum definitions
#include "mint/config.hpp"             // for mint compile-time type definitions
#include "mint/FieldAssociation.hpp"   // for FieldAssociation enum
#include "mint/FieldData.hpp"          // for mint::FieldData
#include "mint/MeshCoordinates.hpp"    // for mint::MeshCoordinates
#include "mint/MeshTypes.hpp"          // for MeshType enum and property traits

#include "slic/slic.hpp"               // for SLIC macros


namespace axom
{

/* Forward declarations */
#ifdef MINT_USE_SIDRE
namespace sidre
{
class Group;
}
#endif

namespace mint
{

/* Forward declarations */
class Field;

/*!
 * \class Mesh
 *
 * \brief Base class that defines the core API common to all Mesh types.
 *
 * A Mesh, \f$ \mathcal{M}(\Omega) \f$, provides an approximation of a physical
 * domain, \f$ \Omega \in \mathcal{R}^d \f$, where \f$ d \in [1,3] \f$ . The
 * Mesh is essentially a discrete representation of a problem and is used to
 * facilitate the analysis, e.g., FEA, CFD, etc. The solution domain is
 * approximated by dividing it into a finite number of <em> nodes </em> and
 * <em> cells </em> at which the variables of the underlying mathematical model
 * (i.e., a PDE) are then computed via a numerical method, such as,
 * Finite Difference (FD), Finite Volume (FV) or Finite Element (FE), chief
 * among them.
 *
 * There are a variety of mesh types. Mint supports the following mesh types:
 *
 *  * <b> Structured Mesh </b> <br />
 *
 *    A <em> structured mesh </em> divides the solution domain according to a
 *    logical grid where each node/cell of the mesh can be uniquely identified
 *    by a corresponding logical ijk-index. The nodes of the mesh are found at
 *    the intersection of the grid lines, but, are explicitly defined via a
 *    mapping onto the physical Cartesian coordinate system of the domain. For
 *    this reason, these types of meshes are also called <em> mapped </em> or
 *    <em> body-fitted </em> meshes.
 *
 *    However, the mesh topology (e.g., connectivity, neighbor information) is
 *    implicitly defined by the logical indexing scheme. For example,  a
 *    structured mesh is composed of <em> quadrilateral </em> cells in 2-D and
 *    <em> hexahedron </em> cells in 3-D. Given the logical index of a cell
 *    one can compute the node indices of the cell and neighbor information by
 *    performing simple shift operations on the associated index.
 *
 *  * <b> Rectilinear Mesh </b> <br />
 *
 *    A <em> rectilinear mesh </em>, also known as a product mesh, is similar
 *    to the <em> structured mesh </em> in that it is also defined according to
 *    a logical indexing scheme and has implicit topology.
 *
 *    However, the nodes and cells on a <em> rectilinear mesh </em> are arranged
 *    on a regular lattice that is axis-aligned with the Cartesian coordinate
 *    system. In contrast to the general <em> structured mesh </em>, the
 *    <em> rectilinear mesh </em> does not explicitly define **all** the nodes
 *    of the mesh. Instead, the nodes are only defined along each coordinate
 *    axis and may have variable spacing between nodes. Given a logical index,
 *    the corresponding physical position of a node can be evaluated by taking
 *    the cartesian product of the corresponding coordinate along each
 *    coordinate axis.
 *
 *  * <b> Uniform Mesh </b> <br />
 *
 *    A <em> uniform mesh </em>, also called a regular mesh, subdivides the
 *    domain in cells that have uniform spacing across each coordinate axis.
 *    Similar to the <em> structured mesh </em>, a <em> uniform mesh </em>
 *    adheres to the same logical indexing scheme and implicit topology
 *    representation. Moreover, The nodes and cells of a <em> uniform </em> mesh
 *    are arranged on a regular lattice as with the <em> rectilinear mesh </em>.
 *    However, both topology and geometry is implicitly defined on a
 *    <em> uniform mesh </em>. The geometry is solely defined by an origin,
 *    \f$ \hat{x_0} \f$ and spacing, \f$ \hat{h} \f$, along each axis. The
 *    coordinates of a node can be evaluated algebraically by the following:
 *    \f$ \hat{p} = \hat{x_0} + \hat{i} \times \hat{h} \f$, where \f$\hat{i}\f$
 *    is the logical ijk-index of the corresponding node.
 *
 *  * <b> Unstructured Mesh </b> <br />
 *
 *    An <em> unstructured mesh </em> stores both node and topology information
 *    explicitly. This allows the flexibility of discretizing the solution
 *    domain using a variety of cell types not just quadrilateral (in 2D) or
 *    hexahedral (in 3D) cells. Due to this added flexibility, the use of
 *    <em> unstructured meshes </em> is more common when dealing with complex
 *    geometries. However, <em> unstructured meshes </em> require additional
 *    storage and generally incur some performance penalty to store, create and
 *    access mesh topology information respectively.
 *
 *    Mint classifies <em> unstructured meshes </em> in two basic types based on
 *    the underlying mesh topology:
 *
 *    * <b> Single Cell Topology </b>
 *
 *      In this case, the <em> unstructured mesh </em> consists of a single cell
 *      type, e.g., a quad or triangle mesh in 2D, or, a hex or tet mesh in 3D.
 *      In this case  the underlying implementation is optimized for the
 *      specified cell type (specified in the constructor).
 *
 *    * <b> Mixed Cell Topology </b>
 *
 *      When <em> mixed cell topology </em> is specified, the <em> unstructured
 *      mesh </em> can be composed of any of the supported cell types, e.g.,
 *      a mesh consisting of both quads and triangles. This mode incurs
 *      additional overhead for storage and access to mesh topology information,
 *      since it requires indirection.
 *
 *    The list of supported cell types for an <em> unstructured mesh </em> is
 *    available in CellTypes.hpp
 *
 *  * <b> Particle Mesh </b> <br />
 *
 *    A <em> particle mesh </em> discretizes the solution domain using
 *    a set of discrete particle elements which correspond to the the nodes
 *    of the mesh. There is no ordering imposed on the particles and the
 *    coordinates of the particles are explicitly defined. A particle mesh
 *    has no connectivity information connecting the particles, which is why
 *    in literature methods using a particle discretization are also referred
 *    to as <em> meshless </em> or <em> meshfree </em> methods.
 *
 * \note Mesh is a base class and cannot be instantiated directly
 *
 * \note Typically, the computational mesh can be defined across one or more
 *  <em> blocks </em>, e.g., for multi-block problems, where each block is
 *  subsequently decomposed into several <em> partitions </em> for parallel
 *  computation. A Mesh instance represents a <em> single </em> partition for
 *  a given block.
 *
 * \see mint::UnstructuredMesh
 * \see mint::StructuredMesh
 * \see mint::CurvilinearMesh
 * \see mint::RectilinearMesh
 * \see mint::UniformMesh
 * \see mint::Field
 * \see mint::FieldData
 * \see mint::MeshTypes
 */
class Mesh
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  Mesh( ) = delete;

  /*!
   * \brief Destructor.
   */
  virtual ~Mesh();

/// \name Mesh Attribute Query Methods
/// @{

  /*!
   * \brief Returns the dimension for this mesh instance.
   * \return ndims the dimension of this mesh instance.
   * \post ndims >= 1 && ndims <= 3
   */
  inline int getDimension() const
  { return m_ndims; }

  /*!
   * \brief Returns the ID of this mesh instance.
   * \return Id the ID of the mesh.
   */
  inline int getBlockId() const
  { return m_block_idx; }

  /*!
   * \brief set the block ID of this mesh instance.
   * 
   * \param [in] ID the new block ID.
   *
   * \post getBlockId() == ID
   */
  void setBlockId( int ID );

  /*!
   * \brief Returns the partition ID of this mesh instance.
   * \return partitionId the partition ID of the mesh.
   */
  inline int getPartitionId() const
  { return m_part_idx; }

  /*!
   * \brief set the partition ID of this mesh instance.
   * 
   * \param [in] ID the new partition ID.
   *
   * \post getPartitionId() == ID
   */
  void setPartitionId( int ID );

  /*!
   * \brief Returns the mesh type of this mesh instance.
   * \return meshType the mesh type
   * \see MeshType
   */
  inline int getMeshType() const
  { return m_type; }

  /*!
   * \brief Returns the number of nodes in this mesh instance.
   * \return N the number of nodes
   * \post N >= 0
   */
  virtual IndexType getNumberOfNodes() const = 0;

  /*!
   * \brief Returns the number of cells in this mesh instance.
   * \return N the number of cells
   * \post N >= 0
   */
  virtual IndexType getNumberOfCells() const = 0;

  /*!
   * \brief Returns the number of faces in this mesh instance.
   * \return N the number of faces
   * \post N >= 0
   */
  virtual IndexType getNumberOfFaces() const = 0;

  /*!
   * \brief Returns the number of edges in this mesh instance.
   * \return N the number of edges
   * \post N >= 0
   */
  virtual IndexType getNumberOfEdges() const = 0;

  /*!
   * \brief Returns the capacity for number of nodes in this mesh instance.
   * \return N the node capacity
   * \post N >= 0
   */
  virtual IndexType getNodeCapacity() const = 0;

  /*!
   * \brief Returns the capacity for number of cell in this mesh instance.
   * \return N the cell capacity
   * \post N >= 0
   */
  virtual IndexType getCellCapacity() const = 0;

  /*!
   * \brief Returns the capacity for number of faces in this mesh instance.
   * \return N the face capacity
   * \post N >= 0
   */
  virtual IndexType getFaceCapacity() const = 0;

  /*!
   * \brief Returns the capacity for number of edges in this mesh instance.
   * \return N the edge capacity
   * \post N >= 0
   */
  virtual IndexType getEdgeCapacity() const = 0;

  /*!
   * \brief Returns the node resize ratio for this mesh instance.
   * \return R the node resize ratio
   */
  virtual double getNodeResizeRatio() const = 0;

  /*!
   * \brief Returns the cell resize ratio for this mesh instance.
   * \return R the cell resize ratio
   */
  virtual double getCellResizeRatio() const = 0;

  /*!
   * \brief Returns the face resize ratio for this mesh instance.
   * \return R the face resize ratio
   */
  virtual double getFaceResizeRatio() const = 0;

  /*!
   * \brief Returns the edge resize ratio for this mesh instance.
   * \return R the edge resize ratio
   */
  virtual double getEdgeResizeRatio() const = 0;

  /*!
   * \brief Checks if this mesh instance has explicit coordinates.
   * \return status true iff the mesh defines coordinates explicitly.
   */
  inline bool hasExplicitCoordinates() const;
  { return m_explicit_coords; }

  /*!
   * \brief Checks if this mesh instance has explicit connectivity.
   * \return status true iff the mesh defines cell connectivity explicitly.
   */
  inline bool hasExplicitConnectivity() const
  { return m_explicit_connectivity; }

  /*!
   * \brief Checks if the mesh has mixed cell types, e.g., consisting of both
   *  triangle and quad elements or hex,pyramid,prisms and tets in 3-D.
   * \return status true iff the mesh has mixed cell types.
   */
  inline bool hasMixedCellTypes() const
  { return m_has_mixed_topology; }

  /*!
   * \brief Checks if this Mesh instance is associated with a Sidre Group.
   * \return status true if the Mesh is associated with a group in a Sidre
   *  hierarchy, else, false.
   */
  inline bool hasSidreGroup() const;

/// @}

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \pre hasExplicitCoordinates()==true
   * \pre dim >= 0 && dim < dimension()
   * \pre dim==X_COORDINATE || dim==Y_COORDINATE || dim==Z_COORDINATE
   * \post ptr != AXOM_NULLPTR
   *
   * \see MeshCoordinates
   */
  /// @{

  virtual double* getCoordinateArray( int dim ) = 0;
  virtual const double* getCoordinateArray( int dim ) const = 0;

  /// @}

  /*!
   * \brief Copy the coordinates of the given node into the provided buffer.
   *
   * \param [in] nodeID the ID of the node in question.
   * \param [in] coords the buffer to copy the coordinates into, of length at
   *  least getDimension().
   *
   * \pre 0 <= nodeID < getNumberOfNodes()
   * \pre coords != AXOM_NULLPTR
   */
  virtual void getNode( IndexType nodeID, double* node ) const = 0;

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored unless hasMixedCellTypes() == true.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getNumberOfCellNodes( IndexType cellID=0 ) const = 0;

  /*!
   * \brief Copy the connectivity of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] cell the buffer into which the connectivity is copied.
   *
   * \return The number of nodes for the given cell.
   * 
   * \pre cell != AXOM_NULLPTR
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType getCell( IndexType cellID, IndexType* cell ) const = 0;

  /*!
   * \brief Return the type of cell this mesh holds. Returns UNDEFINED_CELL if
   *  hasMixedCellTypes() == true.
   */
  virtual CellType getCellType() const = 0;

  /*!
   * \brief Return the type of the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored if hasMixedCellTypes() == false.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual CellType getCellType( IndexType cellID ) const = 0;

/// \name Methods to Create, Access & Remove Fields from a Mesh
/// @{

  /*!
   * \brief Returns const pointer to the FieldData instance with the specified
   *  mesh field association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] association the specified mesh field association
   * \return fd pointer to the requested FieldData instance
   *
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post fd != AXOM_NULLPTR
   *
   * \see FieldAssociation
   * \see FieldData
   */
  inline const FieldData* getFieldData( int association ) const;

  /*!
   * \brief Check if a field with the given name and association exists.
   *
   * \param [in] name the name of the field in query.
   * \param [in] association the field association, e.g., NODE_CENTERED, etc.
   *
   * \return status true if the field exists, else, false.
   *
   * \pre name.empty()==false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool hasField( const std::string& name, int association ) const;

  /*!
   * \brief Creates a new field with the given name and specified mesh field
   *  association, e.g., NODE_CENTERED, CELL_CENTERED, etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] num_components number of components of the field (optional).
   * \param [in] storeInSidre indicates whether to store the field in the
   *  corresponding Sidre group (optional).
   * \param [in] capacity
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post ptr != AXOM_NULLPTR
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         int association,
                         IndexType num_components=1,
                         bool storeInSidre=true );

  /*!
   * \brief Creates a new field from an external buffer that has the given name
   *  and specified mesh field association, e.g., NODE_CENTERED, CELL_CENTERED,
   *  etc.
   *
   * \param [in] name the name of the new field.
   * \param [in] association the mesh field association.
   * \param [in] data pointer to the external data buffer.
   * \param [in] num_components number of components of the field (optional).
   *
   * \return ptr raw pointer to the data buffer of the new field.
   *
   * \pre name.empty() == false
   * \pre hasField( name ) == false
   * \pre data != AXOM_NULLPTR
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   * \post ptr != AXOM_NULLPTR
   * \post ptr == data
   * \post hasField( name ) == true
   *
   * \see FieldAssociation
   */
  template < typename T >
  inline T* createField( const std::string& name,
                         int association,
                         T* data,
                         IndexType num_components=1,
                         IndexType capacity=USE_DEFAULT );

  /*!
   * \brief Removes the field with the given name and specified association.
   *
   * \param [in] name the name of the field to remove.
   * \param [in] association the mesh field association.
   *
   * \return status true if the field is removed successfully, else, false.
   *
   * \pre name.emtpy() == false
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  inline bool removeField( const std::string& name, int association );

  /*!
   * \brief Returns pointer to buffer of the field with the given ane and
   *  specified mesh field association.
   *
   * \param [in] name the name of the requested field.
   * \param [in] association the mesh field association.
   * \param [out] num_components the number of components per tuple (optional).
   *
   * \return ptr raw pointer to the data buffer of the requested field.
   *
   * \pre name.empty() == false
   * \pre hasField( name )
   * \pre association >= 0 && association < NUM_FIELD_ASSOCIATION
   *
   * \see FieldAssociation
   */
  /// @{

  template < typename T >
  inline T* getFieldPtr( const std::string& name,
                         int association,
                         IndexType& num_components );

  template < typename T >
  inline T* getFieldPtr( const std::string& name,
                         int association );

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               int association,
                               IndexType& num_components ) const;

  template < typename T >
  inline const T* getFieldPtr( const std::string& name,
                               int association ) const;

  /// @}

/// @}

/// \name Static Methods
/// @{

#ifdef MINT_USE_SIDRE

  /*!
   * \brief Creates a mesh instance from the given Sidre group.
   *
   * \param [in] group pointer to the root group of the mesh in Sidre.
   * \param [in] topo topology associated with the requested mesh (optional)
   *
   * \return m pointer to a mesh instance corresponding to the specified group.
   *
   * \note If a topology name is not provided, the implementation will construct
   *  a mesh based on the 1st topology group under the parent "topologies"
   *  group.
   *
   * \note Ownership of the resulting mesh object is passed to the caller.
   *
   * \note When using Mint with Sidre, Sidre maintains ownership of the data.
   *  Although the data can be modified via calls to Mint, groups and views
   *  cannot be deleted. The data remains persistent in Sidre once the mesh
   *  object goes out-of-scope.
   *
   * \pre  group != AXOM_NULLPTR
   * \pre  blueprint::validRootGroup( group ) == true
   * \post m != AXOM_NULLPTR
   */
  /// @{
  static Mesh* getMesh( const sidre::Group* group, const std::string& topo );
  static Mesh* getMesh( const sidre::Group* group );
  /// @}

#endif

/// @}

protected:

/// \name Protected Members
/// @{

  int m_ndims;                    /*! mesh dimension */
  int m_type;                     /*! the type of the mesh */
  int m_block_idx;                /*! the Block ID of the mesh */
  int m_part_idx;                 /*! the partition ID of the mesh */

  bool m_explicit_coords;
  bool m_explicit_connectivity;
  bool m_has_mixed_topology;

  FieldData* m_mesh_fields[ NUM_FIELD_ASSOCIATIONS ];

#ifdef MINT_USE_SIDRE
  sidre::Group* m_group;
  std::string m_topology;
#endif

/// @}

/// \name Protected Constructors (used in derived classes )
/// @{

  /*!
   * \brief Mesh Constructor.
   *
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   */
  Mesh( int ndims, int type );

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Constructor for use with a group that already has data.
   *
   * \param [in] group the sidre::Group to use.
   * \param [in] topo optional argument specifying the name of the topology
   *  associated with this Mesh instance.
   *
   * \note If a topology name is not provided, the implementation will construct
   *  a mesh based on the 1st topology group under the parent "topologies"
   *  group.
   *
   * \pre group != AXOM_NULLPTR.
   * \pre blueprint::validRootGroup( group ) == true
   *
   * \see sidre::Group
   */
  /// @{

  Mesh( sidre::Group* group, const std::string& topo );
  explicit Mesh( sidre::Group* group );

  /// @}

  /*!
   * \brief Constructor for use with an empty group.
   *
   * \param [in] ndims the number of dimensions
   * \param [in] type the mesh type.
   * \param [in] blockId the block ID for this mesh instance.
   * \param [in] partId the partition ID for this mesh instance.
   * \param [in] group the sidre::Group to use.
   * \param [in] topo the name of the associated topology group.
   * \param [in] coordset the name of the associated coordset group.
   *
   * \note If a topology and coordset name is not provided a default name is
   *  used by the implementation.
   *
   * \pre group != AXOM_NULLPTR.
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   * \post blueprint::validRootGroup( group )
   *
   * \see sidre::Group
   */
  /// @{

  Mesh( int ndims, int type, int blockId, int partId,
        sidre::Group* group,
        const std::string& topo,
        const std::string& coordset );

  Mesh( int ndims, int type, int blockId, int partId, sidre::Group* group );
  /// @}

  /*!
   * \brief Helper method to return the associated coordset group.
   * \return coordset the associated coordset group.
   *
   * \pre  m_group != AXOM_NULLPTR
   * \pre  blueprint::validRootGroup( m_group )
   * \post blueprint::validCoordsetGroup( coordset )
   */
  sidre::Group* getCoordsetGroup( );

  Mesh( sidre::Group* group, const std::string& topo, const std::string& coordset,
        int ndims, int type );

  /*!
   * \brief Helper method to return the associated topology group.
   * \return topology the associated topology group.
   *
   * \pre  m_group != AXOM_NULLPTR
   * \pre  blueprint::validRootGroup( m_group )
   * \post blueprint::validTopologyGroup( topology )
   */
  sidre::Group* getTopologyGroup();

#endif

/// @}

private:

  /*!
   * \brief Get the info corresponding to the given mesh field association.
   *
   * \param [in] association the mesh field association, e.g., NODE_CENTERED.
   * \param [out] num_tuples the number of tuples in the associated FieldData.
   * \param [out] capacity the capacity of the associated FieldData.
   * \param [out] ratio the resize ratio of the associated FieldData.
   */
  void getFieldInfo( int association, IndexType& num_tuples, 
                     IndexType& capacity, double& ratio ) const;

  /*!
   * \brief Helper method to check if the mesh type is valid.
   * \return status true if the mesh type is valie, else, false.
   */
  inline bool validMeshType( ) const
  { return ( (m_type >= 0) && (mint::NUM_MESH_TYPES) ); }

  /*!
   * \brief Helper method to check if the mesh dimension is valid.
   * \return status true if the mesh dimension is valid, else, false.
   */
  inline bool validDimension( ) const
  { return ( m_ndims >= 1 && m_ndims <= 3 ); }

  /*!
   * \brief Allocates the FieldData internal data-structures.
   * \note Helper method that is called from the constructor.
   */
  void allocateFieldData( );

  /*!
   * \brief Deallocates the FieldData internal data-structures.
   * \note Helper method that is called by the destructor.
   */
  void deallocateFieldData( );

  DISABLE_COPY_AND_ASSIGNMENT( Mesh );
  DISABLE_MOVE_AND_ASSIGNMENT( Mesh );
};

//------------------------------------------------------------------------------
//  IMPLEMENTATION OF TEMPLATE & IN-LINE METHODS
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline bool Mesh::hasSidreGroup( ) const
{
#ifdef MINT_USE_SIDRE
  return ( m_group != AXOM_NULLPTR );
#else
  return false;
#endif
}

//------------------------------------------------------------------------------
inline IndexType Mesh::getNumTuples( int association ) const
{
  switch ( association )
  {
  case NODE_CENTERED:
    return getNumberOfNodes();
  case CELL_CENTERED:
    return getNumberOfCells();
  case FACE_CENTERED:
    return getNumberOfFaces();
  default:
    SLIC_ASSERT( association==EDGE_CENTERED );
    return getNumberOfEdges();
  } // END switch
}

//------------------------------------------------------------------------------
inline IndexType Mesh::getCapacity( int association ) const
{
  switch ( association )
  {
  case NODE_CENTERED:
    return getNodeCapacity();
  case CELL_CENTERED:
    return getCellCapacity();
  case FACE_CENTERED:
    return getFaceCapacity();
  default:
    SLIC_ASSERT( association == EDGE_CENTERED );
    return getEdgeCapacity();
  } // END switch
}

//------------------------------------------------------------------------------
inline double Mesh::getResizeRatio( int association ) const
{
  switch ( association )
  {
  case NODE_CENTERED:
    return getNodeResizeRatio();
  case CELL_CENTERED:
    return getCellResizeRatio();
  case FACE_CENTERED:
    return getFaceResizeRatio();
  default:
    SLIC_ASSERT( association == EDGE_CENTERED );
    return getEdgeResizeRatio();
  } // END switch
}

//------------------------------------------------------------------------------
inline const FieldData* Mesh::getFieldData( int association ) const
{
  SLIC_ERROR_IF( association < 0 || association >= NUM_FIELD_ASSOCIATIONS,
                  "invalid field association [" << association << "]" );
  SLIC_ERROR_IF( m_mesh_fields[ association ]==AXOM_NULLPTR,
              "null field data object w/association [" << association << "]" );
  SLIC_ERROR_IF( m_type==PARTICLE_MESH && association != NODE_CENTERED,
              "a particle mesh may only store node-centered fields" );

  return m_mesh_fields[ association ];
}

//------------------------------------------------------------------------------
inline bool Mesh::hasField( const std::string& name, int association ) const
{
  const FieldData* fd = getFieldData( association );
  SLIC_ASSERT( fd != AXOM_NULLPTR );
  return fd->hasField( name );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::createField( const std::string& name,
                             int association,
                             IndexType num_components,
                             bool storeInSidre )
{
  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = getNumTuples( association );
  IndexType capacity = getCapacity( association );
  double ratio = getResizeRatio( association );
  T* ptr = fd->createField< T >( name, num_tuples, num_components, capacity, 
                                                          storeInSidre, ratio );
  if ( num_tuples > 0 ) 
  {
    SLIC_ASSERT( ptr != AXOM_NULLPTR );
  }

  return ( ptr );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::createField( const std::string& name,
                             int association,
                             T* data,
                             IndexType num_components,
                             IndexType capacity )
{
  SLIC_ASSERT( data != AXOM_NULLPTR );

  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = getNumTuples( association );
  T* ptr = fd->createField< T >( name, data, num_tuples, num_components, capacity );
  SLIC_ASSERT( ptr == data );

  return ( ptr );
}

//------------------------------------------------------------------------------
inline bool Mesh::removeField( const std::string& name, int association )
{
  bool status   = false;
  FieldData* fd = const_cast< FieldData* >( getFieldData( association ) );

  const bool hasField = fd->hasField( name );
  SLIC_WARNING_IF( !hasField, "field [" << name << "] does not exist!" );

  if ( hasField )
  {
    fd->removeField( name );
    status = true;
  }

  return ( status );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::getFieldPtr( const std::string& name,int association )
{
  IndexType num_components = 0;
  return getFieldPtr< T >( name, association, num_components );
}

//------------------------------------------------------------------------------
template < typename T >
inline T* Mesh::getFieldPtr( const std::string& name,
                             int association,
                             IndexType& num_components )
{
  const Mesh* self = const_cast< const Mesh* >( this );
  const T* ptr = self->getFieldPtr< T >( name, association, num_components );
  return ( const_cast< T* >( ptr ) );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* Mesh::getFieldPtr( const std::string& name,
                                   int association ) const
{
  IndexType num_components = 0;
  return getFieldPtr< T >( name, association, num_components );
}

//------------------------------------------------------------------------------
template < typename T >
inline const T* Mesh::getFieldPtr( const std::string& name,
                             int association,
                             IndexType& num_components ) const
{
  const FieldData* fd = getFieldData( association );
  SLIC_ASSERT( fd != AXOM_NULLPTR );

  IndexType num_tuples = 0;
  const T* ptr = fd->getFieldPtr< T >( name, num_tuples, num_components );
  SLIC_ASSERT( ptr != AXOM_NULLPTR );

  return ( ptr );
}

} /* namespace mint */
} /* namespace axom */

#endif /* MESH_HXX_ */
