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

#ifndef MINT_STRUCTUREDMESH_HPP_
#define MINT_STRUCTUREDMESH_HPP_

#include "axom/core/Types.hpp"            // for axom types
#include "axom/core/Macros.hpp"           // for axom macros

#include "axom/mint/config.hpp"           // for compile-time definitions
#include "axom/mint/mesh/CellTypes.hpp"   // for the CellTypes enum
#include "axom/mint/mesh/Mesh.hpp"        // for mint::Mesh base class

#include "axom/slic/interface/slic.hpp"   // for SLIC macros

#include <cstring>                        // for std::memset

namespace axom
{
namespace mint
{

/*!
 * \class StructuredMesh
 *
 * \brief Base class that defines the core API common for structured mesh types.
 *
 *  The StructuredMesh class derives from the abstract Mesh base class and
 *  implements the core API for Structured meshes. Specifically, the
 *  StrucrturedMesh  defines and implements all the extent-based operations.
 *
 * \see UniformMesh
 * \see RectilinearMesh
 * \see CurvilinearMesh
 * \see Mesh
 */
class StructuredMesh : public Mesh
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  StructuredMesh() = delete;

/// \name Virtual methods
/// @{

  /*!
   * \brief Destructor.
   */
  virtual ~StructuredMesh()
  {}

/// \name Cells
/// @{

  /*!
   * \brief Return the number of cells in the mesh.
   */
  virtual IndexType getNumberOfCells() const final override
  {
    IndexType n_cells = 1;
    for ( int dim = 0; dim < m_ndims; ++dim )
    {
      n_cells *= getCellExtent( dim );
    }

    return n_cells;
  }

  /*!
   * \brief Return the type of cell this mesh holds. SEGMENT, QUAD, or HEX
   *  depending on the dimension.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   */
  virtual CellType
  getCellType( IndexType AXOM_NOT_USED(cellID)=0 ) const final override
  {
    return ( m_ndims == 1 ) ? SEGMENT :
           ( m_ndims == 2 ) ? QUAD : HEX;
  }

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType
  getNumberOfCellNodes( IndexType AXOM_NOT_USED(cellID)=0 ) const final override
  {
    return ( m_ndims == 1 ) ? 2 :
           ( m_ndims == 2 ) ? 4 : 8;
  }

  /*!
   * \brief Copy the connectivity of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] nodes the buffer into which the connectivity is copied, must
   *  be of length at least getNumberOfCellNodes().
   *
   * \return The number of nodes for the given cell.
   *
   * \pre nodes != AXOM_NULLPTR
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType
  getCellNodes( IndexType cellID, IndexType* nodes ) const final override;

  virtual IndexType
  getNumberOfCellFaces( IndexType AXOM_NOT_USED(cellID)=0 ) const final override
  {
    CellType cell_type = getCellType();
    return getCellInfo( cell_type ).num_faces;
  }

  virtual IndexType
  getCellFaces( IndexType cellID, IndexType* faces ) const final override;

/// @}

/// \name Nodes
/// @{

  /*!
   * \brief Return the number of nodes in the mesh.
   */
  virtual IndexType getNumberOfNodes() const final override
  {
    IndexType n_nodes = 1;
    for ( int dim = 0; dim < m_ndims; ++dim )
    {
      n_nodes *= getNodeExtent( dim );
    }

    return n_nodes;
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
  virtual void getNode( IndexType nodeID, double* node ) const override = 0;

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \note if hasExplicitCoordinates() == true then the length of the returned
   *  buffer is getNumberOfNodes(). Otherwise the UniformMesh returns
   *  nullptr and the RectilinearMesh returns a pointer to the associated
   *  dimension scale which is of length
   *  static_cast< RectilinearMesh* >( this )->getNodeExtent( dim ).
   *
   * \pre dim >= 0 && dim < dimension()
   * \pre dim == X_COORDINATE || dim == Y_COORDINATE || dim == Z_COORDINATE
   */
  /// @{

  virtual double* getCoordinateArray( int dim ) override = 0;
  virtual const double* getCoordinateArray( int dim ) const override = 0;

  /// @}

/// @}

/// \name Faces
/// @{

  /*!
   * \brief Return the number of faces in the mesh.
   */
  virtual IndexType getNumberOfFaces() const final override
  { return getTotalNumFaces( 0 ) + getTotalNumFaces( 1 ) + getTotalNumFaces( 2 ); }

  virtual CellType
  getFaceType( IndexType AXOM_NOT_USED(cellID)=0 ) const final override
  {
    return ( m_ndims == 2 ) ? SEGMENT :
           ( m_ndims == 3 ) ? QUAD : UNDEFINED_CELL;
  }

  virtual IndexType 
  getNumberOfFaceNodes( IndexType AXOM_NOT_USED(faceID)=0 ) const final override
  {
    return ( m_ndims == 2 ) ? 2 :
           ( m_ndims == 3 ) ? 4 : 0;
  }

  virtual IndexType
  getFaceNodes( IndexType faceID, IndexType* nodes ) const final override;

  virtual void 
  getFaceCells( IndexType faceID, IndexType& cellIDOne, 
                IndexType& cellIDTwo ) const final override;

/// @}

/// \name Edges
/// @{

  /*!
   * \brief Return the number of edges in the mesh.
   */
  virtual IndexType getNumberOfEdges() const final override
  { return m_num_edges; }

/// @}

  /*!
   * \brief Returns true iff the mesh was constructed with external arrays.
   * \return status true if the mesh points to external buffers, else, false.
   */
  virtual bool isExternal() const override
  { return false; }

/// @}

/// \name Data Accessor Methods
/// @{

/// \name Cells
/// @{

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre getDimension() == 2.
   */
  inline void getCellNodes( IndexType i, IndexType j, IndexType* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j,k)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre getDimension() == 3.
   */
  inline void
  getCellNodes( IndexType i, IndexType j, IndexType k, IndexType* cell) const;

  inline void getCellFaces( IndexType i, IndexType j, IndexType* faces ) const
  {
    const IndexType cellID = getCellLinearIndex( i, j );
    getCellFacesInternal( cellID, j, faces );
  }

  inline void 
  getCellFaces( IndexType i, IndexType j, IndexType k, IndexType* faces ) const
  {
    const IndexType cellID = getCellLinearIndex( i, j, k );
    getCellFacesInternal( cellID, j, k, faces );
  }

/// @}

/// \name Faces
/// @{

  inline IndexType getIFaceNodes( IndexType faceID, IndexType* nodes ) const;

  inline IndexType getJFaceNodes( IndexType faceID, IndexType* nodes ) const;

  inline IndexType getKFaceNodes( IndexType faceID, IndexType* nodes ) const;

  inline void
  getIFaceCells( IndexType faceID, IndexType& cellIDOne, 
                 IndexType& cellIDTwo ) const;
  inline void
  getJFaceCells( IndexType faceID, IndexType& cellIDOne, 
                 IndexType& cellIDTwo ) const;
  inline void
  getKFaceCells( IndexType faceID, IndexType& cellIDOne, 
                 IndexType& cellIDTwo ) const;

/// @}

/// @}


/// \name Attribute Querying Methods
/// @{

  /*!
   * \brief Returns the number of nodes along the given dimension.
   * \param [in] dim the dimension to querry.
   * \pre 0 <= dim < 3
   */
  inline IndexType getNodeExtent( IndexType dim ) const
  {
    SLIC_ASSERT( 0 <= dim && dim < 3 );
    return m_node_extent[ dim ];
  }

  inline 
  void getGlobalNodeExtent( IndexType dim, int64& low, int64& high ) const
  {
    SLIC_ASSERT( 0 <= dim && dim < 3 );
    low = m_global_node_extent[ 2 * dim ];
    high = m_global_node_extent[ 2 * dim + 1 ];
  }

  
  void setGlobalNodeExtent( int ndims, const int64* extent );
  
  /*!
   * \brief Returns the number of nodes along the given dimension.
   * \param [in] dim the dimension to querry.
   * \pre 0 <= dim < 3
   */
  inline IndexType getCellExtent( IndexType dim ) const
  {
    SLIC_ASSERT( 0 <= dim && dim < 3 );
    return m_cell_extent[ dim ];
  }

  inline IndexType getTotalNumFaces( int direction ) const
  {
    SLIC_ASSERT( I_DIRECTION <= direction );
    SLIC_ASSERT( direction <= K_DIRECTION );
    return m_total_faces[ direction ];
  }

/// @}

/// \name Indexing Helper Methods
/// @{

/// \name Nodes
/// @{

  /*!
   * \brief Returns stride to the second dimension of nodes.
   * \post jp >= 0.
   */
  inline IndexType nodeJp() const
  { return getNodeExtent( I_DIRECTION ); }

  /*!
   * \return kp stride to the third dimension of nodes.
   * \post kp >= 0.
   */
  inline IndexType nodeKp() const
  { return m_node_kp; }

  /*!
   * \brief Returns the linear index corresponding to the given logical node
   *  indices.
   * \param [in] i logical node index of the first dimension.
   * \param [in] j logical node index of the second dimension.
   * \param [in] k logical node index of the third dimension (optional)
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   */
  /// @{
  inline IndexType getNodeLinearIndex( IndexType i, IndexType j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    return i + j * nodeJp();
  }

  inline IndexType 
  getNodeLinearIndex( IndexType i, IndexType j, IndexType k  ) const
  { 
    SLIC_ASSERT( m_ndims == 3 );
    return i + j * nodeJp() + k * m_node_kp;
  }
  /// @}

  /*!
   * \brief Given the 1D linear index of a node, this method computes the
   *  corresponding i-j-k grid index.
   *
   * \param [in] linearIdx the local flat index.
   * \param [out] i the corresponding grid index along the I_DIRECTION.
   * \param [out] j the corresponding grid index along the J_DIRECTION.
   * \param [out] k the corresponding grid index along the K_DIRECTION.
   *
   * \note The computed i-j-k grid indices are expected to be the shifted
   *  topological coordinates within the local frame of reference for the
   *  extent.
   *
   * \pre linearIdx >= 0 && linearIdx < getNumNodes()
   * \post i >= 0 && i < size( I_DIRECTION )
   * \post j >= 0 && j < size( J_DIRECTION )
   * \post k >= 0 && k < size( K_DIRECTION )
   */
  /// @{
  inline void 
  getNodeGridIndex( IndexType linearIdx, IndexType& i, IndexType& j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    j = linearIdx / nodeJp();
    i = linearIdx - j * nodeJp();
  }

  inline void 
  getNodeGridIndex( IndexType linearIdx, IndexType& i, IndexType& j, 
                    IndexType& k ) const
  {
    SLIC_ASSERT( m_ndims == 3 );
    k = linearIdx / m_node_kp;
    const IndexType temp = linearIdx - k * m_node_kp;
    j = temp / nodeJp();
    i = temp - j * nodeJp();
  }
  /// @}

/// @}

/// \name Cells
/// @{

  /*!
   * \brief Returns stride to the second dimension of cells.
   * \post jp >= 0.
   */
  inline IndexType cellJp() const
  { return getCellExtent( I_DIRECTION ); }

  /*!
   * \brief Returns stride to the third dimension of cells.
   * \post kp >= 0.
   */
  inline IndexType cellKp() const
  { return m_cell_kp; }

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   */
  /// @{
  inline IndexType getCellLinearIndex( IndexType i, IndexType j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    return i + j * cellJp();
  }

  inline IndexType
  getCellLinearIndex( IndexType i, IndexType j, IndexType k ) const
  {
    SLIC_ASSERT( m_ndims == 3 );
    return i + j * cellJp()  + k * m_cell_kp;
  }
  /// @}

  /*!
   * \brief Given the 1D linear index of a cell, this method computes
   *  the corresponding i-j-k grid index.
   *
   * \param [in] linearIdx the local flag index.
   * \param [out] i the corresponding grid index along the I_DIRECTION.
   * \param [out] j the corresponding grid index along the J_DIRECTION.
   * \param [out] k the corresponding grid index along the K_DIRECTION.
   *
   * \note The computed i-j-k grid indices are expected to be the shifted
   *  topological coordinates within the local frame of reference for the
   *  extent.
   *
   * \pre linearIdx >= 0 && linearIdx < getNumNodes()
   * \post i >= 0 && i < size( I_DIRECTION )-1
   * \post j >= 0 && j < size( J_DIRECTION )-1
   * \post k >= 0 && k < size( K_DIRECTION )-1
   */
  /// @{
  inline void 
  getCellGridIndex( IndexType linearIdx, IndexType& i, IndexType& j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    i = linearIdx % cellJp();
    j = linearIdx / cellJp();
  }

  inline void 
  getCellGridIndex( IndexType linearIdx, IndexType& i, IndexType& j,
                    IndexType& k ) const
  {
    SLIC_ASSERT( m_ndims == 3 );
    k = linearIdx / m_cell_kp;
    const IndexType temp = linearIdx - k * m_cell_kp;
    j = temp / cellJp();
    i = temp - j * cellJp();
  }
  /// @}

/// @}

/// @}

protected:

  /*!
   * \brief Constructs a structured mesh instance from the given extent.
   * \param [in] meshType the mesh type
   * \param [in] dimension the mesh dimension
   * \param [in] cell_ext the structured mesh node extent.
   */
  StructuredMesh( int meshType, int dimension, const IndexType* node_ext );

  StructuredMesh( int meshType, IndexType Ni, IndexType Nj, IndexType Nk );

#ifdef MINT_USE_SIDRE

  /*!
   * \brief Constructs a structured mesh instance from the given group.
   *
   * \param [in] group pointer to the sidre::Group
   * \param [in] topo the topology name to use, an empty string may be supplied.
   *
   * \note If an empty string is supplied for the topology name, the code will
   *  use the 1st topology group under the parent "topologies" group.
   *
   * \pre group !=  nullptr
   * \pre blueprint::isValidRootGroup( group ) == true
   *
   * \note This constructor forwards this call to the parent Mesh class.
   *
   * \see Mesh( sidre::Group* group, const std::string& topo )
   */

  StructuredMesh( sidre::Group* group, const std::string& topo );

  /*!
   * \brief Constructs a structured mesh instance on the specified group.
   *
   * \param [in] meshType the mesh type
   * \param [in] dimension the dimension of the mesh
   * \param [in] group pointer to the group in the Sidre hierarchy.
   * \param [in] topo the topology name to use, may be an empty string.
   * \param [in] coordset the coordset name to use, may be an empty string.
   *
   * \note If an empty string is supplied for the topology and coordset name
   *  respectively, an internal default name will be provided by the
   *  implementation.
   *
   * \pre 1 <= dimension <= 3
   * \pre group != nullptr
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   *
   * \post blueprint::isValidRootGroup( group )
   *
   * \note This constructor forwards this call to the parent Mesh class.
   *
   * \see Mesh( int ndims, int type, sidre::Group*,
   *            const std::string& topo, const std::string& coordset );
   */
  StructuredMesh( int meshType, int dimension, const IndexType* node_ext,
                  sidre::Group* group,
                  const std::string& topo,
                  const std::string& coordset );

  StructuredMesh( int meshType, IndexType Ni, IndexType Nj, IndexType Nk, 
                  sidre::Group* group, const std::string& topo,
                  const std::string& coordset );

#endif

  inline IndexType
  getCellFacesInternal( IndexType cellID, IndexType j, IndexType* faces ) const;

  inline IndexType
  getCellFacesInternal( IndexType cellID, IndexType j, IndexType k,
                        IndexType* faces ) const;

  inline IndexType getIFaceKIndex( IndexType faceID ) const
  { return faceID / (getNodeExtent( 0 ) * getCellExtent( 1 )); }

  inline void getIFaceGridIndex( IndexType faceID, IndexType& i, IndexType& j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    SLIC_ASSERT( 0 <= faceID < getTotalNumFaces( 0 ) );

    j = faceID / getNodeExtent( 0 );
    i = faceID - getNodeExtent( 0 ) * j;
  }

  inline void getIFaceGridIndex( IndexType faceID, IndexType& i, IndexType& j, 
                                 IndexType& k ) const
  {
    SLIC_ASSERT( m_ndims == 3 );
    SLIC_ASSERT( 0 <= faceID < getTotalNumFaces( 0 ) );

    k = getIFaceKIndex( faceID );
    const IndexType temp = (getNodeExtent( 0 ) * getCellExtent( 1 )) * k;
    j = (faceID - temp) / getNodeExtent( 0 );
    i = faceID - getNodeExtent( 0 ) * j - temp;
  }

  inline IndexType shiftJFaceID( IndexType faceID ) const
  { return faceID - getTotalNumFaces( 0 ); }

  inline IndexType getJFaceKIndex( IndexType shiftedID ) const
  { return shiftedID / (getCellExtent( 0 ) * getNodeExtent( 1 )); }

  inline IndexType getJFaceJIndex( IndexType shiftedID, IndexType k ) const
  {
    const IndexType k_stride = getCellExtent( 0 ) * getNodeExtent( 1 );
    return (shiftedID - k_stride * k) / getCellExtent( 0 );
  }

  inline void getJFaceGridIndex( IndexType shiftedID, IndexType& i, IndexType& j ) const
  {
    SLIC_ASSERT( m_ndims == 2 );
    SLIC_ASSERT( 0 <= shiftedID && shiftedID < getTotalNumFaces( 1 ) );

    j = shiftedID / getCellExtent( 0 );
    i = shiftedID - getCellExtent( 0 ) * j;
  }

  inline void getJFaceGridIndex( IndexType shiftedID, IndexType& i, IndexType& j, 
                                 IndexType& k ) const 
  {
    SLIC_ASSERT( m_ndims == 3 );
    SLIC_ASSERT( 0 <= shiftedID && shiftedID < getTotalNumFaces( 1 ) );

    k = getKFaceKIndex( shiftedID );
    j = getKFaceJIndex( shiftedID, k );
    i = shiftedID - getCellExtent( 0 ) * j - (getCellExtent( 0 ) * getNodeExtent ( 1 )) * k;
  }
  inline IndexType shiftKFaceID( IndexType faceID ) const
  { return faceID - getTotalNumFaces( 0 ) - getTotalNumFaces( 1 ); }

  inline IndexType getKFaceKIndex( IndexType shiftedID ) const
  { return shiftedID / cellKp(); }

  inline IndexType getKFaceJIndex( IndexType shiftedID, IndexType k ) const
  { return (shiftedID - cellKp() * k) / cellJp(); }

  inline void getKFaceGridIndex( IndexType shiftedID, IndexType& i, IndexType& j, 
                                 IndexType& k ) const 
  {
    SLIC_ASSERT( m_ndims == 3 );
    SLIC_ASSERT( 0 <= shiftedID && shiftedID < getTotalNumFaces( 2 ) );

    k = getKFaceKIndex( shiftedID );
    j = getKFaceJIndex( shiftedID, k );
    i = shiftedID - cellJp() * j - cellKp() * k;
  }

  void structuredInit();

  IndexType m_node_extent[ 3 ] = { 0, 0, 0 };
  int64 m_global_node_extent[ 6 ] = { 0, 0, 0, 0, 0, 0 };

  IndexType m_node_jp = 0;
  IndexType m_node_kp = 0;

  IndexType m_cell_extent[ 3 ] = { 0, 0, 0 };
  IndexType m_cell_jp = 0;
  IndexType m_cell_kp = 0;
  IndexType m_cell_node_offsets[ 8 ] = { 0, 0, 0, 0, 0, 0, 0, 0 };

  IndexType m_total_faces[ 3 ] = { 0, 0, 0 };

  IndexType m_num_edges = 0;

private:
  DISABLE_COPY_AND_ASSIGNMENT( StructuredMesh );
  DISABLE_MOVE_AND_ASSIGNMENT( StructuredMesh );
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline
IndexType StructuredMesh::getCellNodes( IndexType cellID, IndexType* nodes ) const
{
  SLIC_ASSERT( nodes != AXOM_NULLPTR );
  SLIC_ASSERT( 0 <= cellID && cellID < getNumberOfCells() );

  const IndexType num_cell_nodes = getNumberOfCellNodes();

  // Calculate logical indices of the cell's first corner node.
  IndexType ii, jj, kk;
  getCellGridIndex( cellID, ii, jj, kk );

  // Use the offsets table to get the all the cell nodes.
  const IndexType n0 = getNodeLinearIndex( ii, jj, kk );
  for ( IndexType i = 0 ; i < num_cell_nodes ; ++i )
  {
    nodes[ i ] = n0 + m_cell_node_offsets[ i ];
  }
  return num_cell_nodes;
}

//------------------------------------------------------------------------------
inline IndexType
StructuredMesh::getCellFaces( IndexType cellID, IndexType* faces ) const
{
  SLIC_ASSERT( faces != AXOM_NULLPTR );
  SLIC_ASSERT( 0 <= cellID && cellID < getNumberOfCells() );

  if ( m_ndims == 2 )
  {
    IndexType j = cellID / cellJp();
    return getCellFacesInternal( cellID, j, faces );
  }
  else if ( m_ndims == 3 )
  {
    IndexType k = cellID / m_cell_kp;
    IndexType j = (cellID - k * m_cell_kp) / cellJp();
    return getCellFacesInternal( cellID, k, j, faces );
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------
inline IndexType
StructuredMesh::getFaceNodes( IndexType faceID, IndexType* nodes ) const
{
  SLIC_ASSERT( 0 <= faceID && faceID < getNumberOfFaces() );

  /* Check which kind of face is being asked for */
  if ( faceID < getTotalNumFaces( 0 ) )
  {
    /* It's a I_DIRECTION face */
    return getIFaceNodes( faceID, nodes );
  }
  else if ( faceID < getTotalNumFaces( 0 ) + getTotalNumFaces( 2 ) )
  {
    /* It's a J_DIRECTION face */
    return getJFaceNodes( faceID, nodes );
  }
  else
  {
    /* It's a K_DIRECTION face */
    return getKFaceNodes( faceID, nodes );
  }
}

//------------------------------------------------------------------------------
inline void
StructuredMesh::getFaceCells( IndexType faceID, IndexType& cellIDOne, 
                              IndexType& cellIDTwo ) const
{
  /* Check which kind of face is being asked for */
  if ( faceID < getTotalNumFaces( 0 ) )
  {
    /* It's a I_DIRECTION face */
    return getIFaceCells( faceID, cellIDOne, cellIDTwo );
  }
  else if ( faceID < getTotalNumFaces( 0 ) + getTotalNumFaces( 2 ) )
  {
    /* It's a J_DIRECTION face */
    return getJFaceCells( faceID, cellIDOne, cellIDTwo );
  }
  else
  {
    /* It's a K_DIRECTION face */
    return getKFaceCells( faceID, cellIDOne, cellIDTwo );
  }
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCellNodes( IndexType i, IndexType j,
                                          IndexType* nodes ) const
{
  SLIC_ASSERT( getDimension() == 2 );
  SLIC_ASSERT( nodes != AXOM_NULLPTR );

  /* Get the node index of the first node in the cell. */
  const IndexType n0 = getNodeLinearIndex(i, j);
  
  /* Fill in the node connectivity with the known offsets */
  for ( IndexType ii = 0 ; ii < 4 ; ++ii )
  {
    nodes[ ii ] = n0 + m_cell_node_offsets[ ii ];
  }
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCellNodes( IndexType i, IndexType j, IndexType k,
                                          IndexType* nodes ) const
{
  SLIC_ASSERT( getDimension() == 3 );
  SLIC_ASSERT( nodes != AXOM_NULLPTR );

  /* Get the node index of the first node in the cell. */
  const IndexType n0 = getNodeLinearIndex(i, j, k);

  /* Fill in the node connectivity with the known offsets */
  for ( IndexType ii = 0 ; ii < 8 ; ++ii )
  {
    nodes[ ii ] = n0 + m_cell_node_offsets[ ii ];
  }
}

//------------------------------------------------------------------------------
inline IndexType 
StructuredMesh::getIFaceNodes( IndexType faceID, IndexType* nodes ) const
{
  SLIC_ASSERT( 0 <= faceID && faceID < getTotalNumFaces( 0 ) );

  if ( m_ndims == 2 )
  {
    nodes[ 0 ] = faceID;
    nodes[ 1 ] = nodes[ 0 ] + m_cell_node_offsets[ 3 ];
    return 2;
  }
  else if ( m_ndims == 3 )
  {
    SLIC_ASSERT( m_ndims == 3 );
    const IndexType k = getIFaceKIndex( faceID );
    nodes[ 0 ] = faceID + getCellExtent( 0 ) * k;
    nodes[ 1 ] = nodes[ 0 ] + m_cell_node_offsets[ 3 ];
    nodes[ 2 ] = nodes[ 0 ] + m_cell_node_offsets[ 7 ];
    nodes[ 3 ] = nodes[ 0 ] + m_cell_node_offsets[ 4 ];
    return 4;
  }
  else
  {
    return 0;
  }
}

//------------------------------------------------------------------------------
inline IndexType 
StructuredMesh::getJFaceNodes( IndexType faceID, IndexType* nodes ) const
{
  SLIC_ASSERT( getTotalNumFaces( 0 ) <= faceID && 
               faceID < getTotalNumFaces( 0 ) + getTotalNumFaces( 1 ) );

  const IndexType shiftedID = shiftJFaceID( faceID );
  if ( m_ndims == 2 )
  {
    const IndexType j = getJFaceJIndex( shiftedID, 0 );
    nodes[ 0 ] = shiftedID + j;
    nodes[ 1 ] = nodes[ 0 ] + 1;
    return 2;
  }
  else
  {
    SLIC_ASSERT( m_ndims == 3 );
    const IndexType k = getJFaceKIndex( shiftedID );
    const IndexType j = getJFaceJIndex( shiftedID, k );
    nodes[ 0 ] = shiftedID + j + getNodeExtent( 1 ) * k;
    nodes[ 1 ] = nodes[ 0 ] + 1;
    nodes[ 2 ] = nodes[ 0 ] + m_cell_node_offsets[ 5 ];
    nodes[ 3 ] = nodes[ 0 ] + m_cell_node_offsets[ 4 ];
    return 4;
  }
}

//------------------------------------------------------------------------------
inline IndexType 
StructuredMesh::getKFaceNodes( IndexType faceID, IndexType* nodes ) const
{
  SLIC_ASSERT( getTotalNumFaces( 0 ) + getTotalNumFaces( 1 ) <= faceID && 
               faceID < getNumberOfFaces() );
  SLIC_ASSERT( m_ndims == 3 );

  const IndexType shiftedID = shiftKFaceID( faceID );
  const IndexType k = getKFaceKIndex( shiftedID );
  const IndexType j = getKFaceJIndex( shiftedID, k );
  nodes[ 0 ] = shiftedID + j + (getCellExtent( 0 ) + getCellExtent( 1 ) + 1) * k;
  nodes[ 1 ] = nodes[ 0 ] + 1;
  nodes[ 2 ] = nodes[ 0 ] + m_cell_node_offsets[ 5 ];
  nodes[ 3 ] = nodes[ 0 ] + m_cell_node_offsets[ 4 ];

  return 4;
}

//------------------------------------------------------------------------------
inline void
StructuredMesh::getIFaceCells( IndexType faceID, IndexType& cellIDOne, 
                               IndexType& cellIDTwo ) const
{
  SLIC_ASSERT( m_ndims == 2 || m_ndims == 3 );

  IndexType i, j;
  if ( m_ndims == 2 )
  {
    getIFaceGridIndex( faceID, i, j );
    
    if ( i == 0 )
    {
      /* The first cell doesn't exist */
      cellIDOne = -1;
      cellIDTwo = getCellLinearIndex( i, j );
    }
    else if ( i == getNodeExtent( 0 ) )
    {
      /* The second cell doesn't exist */
      cellIDOne = getCellLinearIndex( i - 1, j );
      cellIDTwo = -1;
    }
    else
    {
      /* Both cells exist */
      cellIDOne = getCellLinearIndex( i - 1, j );
      cellIDTwo = cellIDOne + 1;
    }
  }
  else
  {
    IndexType k;
    getIFaceGridIndex( faceID, i, j, k );
    
    if ( i == 0 )
    {
      /* The first cell doesn't exist */
      cellIDOne = -1;
      cellIDTwo = getCellLinearIndex( i, j, k );
    }
    else if ( i == getNodeExtent( 0 ) )
    {
      /* The second cell doesn't exist */
      cellIDOne = getCellLinearIndex( i - 1, j, k );
      cellIDTwo = -1;
    }
    else
    {
      /* Both cells exist */
      cellIDOne = getCellLinearIndex( i - 1, j, k );
      cellIDTwo = cellIDOne + 1;
    }
  }
}

//------------------------------------------------------------------------------
inline void
StructuredMesh::getJFaceCells( IndexType faceID, IndexType& cellIDOne, 
                               IndexType& cellIDTwo ) const
{
  SLIC_ASSERT( m_ndims == 2 || m_ndims == 3 );

  const IndexType shiftedID = shiftJFaceID( faceID );
  IndexType i , j;
  if ( m_ndims == 2 )
  {
    getJFaceGridIndex( shiftedID, i, j );
    
    if ( j == 0 )
    {
      /* The first cell doesn't exist */
      cellIDOne = -1;
      cellIDTwo = getCellLinearIndex( i, j );
    }
    else if ( j == getNodeExtent( 1 ) )
    {
      /* The second cell doesn't exist */
      cellIDOne = getCellLinearIndex( i , j - 1 );
      cellIDTwo = -1;
    }
    else
    {
      /* Both cells exist */
      cellIDOne = getCellLinearIndex( i , j - 1 );
      cellIDTwo = cellIDOne + cellJp();
    }
  }
  else
  {
    IndexType k;
    getJFaceGridIndex( shiftedID, i, j , k );
    
    if ( j == 0 )
    {
      /* The first cell doesn't exist */
      cellIDOne = -1;
      cellIDTwo = getCellLinearIndex( i, j, k );
    }
    else if ( j == getNodeExtent( 1 ) )
    {
      /* The second cell doesn't exist */
      cellIDOne = getCellLinearIndex( i, j - 1, k );
      cellIDTwo = -1;
    }
    else
    {
      /* Both cells exist */
      cellIDOne = getCellLinearIndex( i, j - 1, k );
      cellIDTwo = cellIDOne + cellJp();
    }
  }
}

//------------------------------------------------------------------------------
inline void
StructuredMesh::getKFaceCells( IndexType faceID, IndexType& cellIDOne, 
                               IndexType& cellIDTwo ) const
{
  SLIC_ASSERT( m_ndims == 3 );

  const IndexType shiftedID = shiftKFaceID( faceID );
  IndexType i , j, k;
  getJFaceGridIndex( shiftedID, i, j , k );
  
  if ( k == 0 )
  {
    /* The first cell doesn't exist */
    cellIDOne = -1;
    cellIDTwo = getCellLinearIndex( i, j, k );
  }
  else if ( k == getNodeExtent( 2 ) )
  {
    /* The second cell doesn't exist */
    cellIDOne = getCellLinearIndex( i, j, k - 1 );
    cellIDTwo = -1;
  }
  else
  {
    /* Both cells exist */
    cellIDOne = getCellLinearIndex( i, j, k - 1 );
    cellIDTwo = cellIDOne + cellKp();
  }
}


//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellFacesInternal( IndexType cellID, 
                                                       IndexType j, 
                                                       IndexType* faces ) const
{
  SLIC_ASSERT( m_ndims == 2 );
  SLIC_ASSERT( 0 <= cellID && cellID < getNumberOfCells() );
  SLIC_ASSERT( 0 <= j && j < getCellExtent( 1 ) );
  SLIC_ASSERT( faces != AXOM_NULLPTR );

  /* The I_DIRECTION faces */
  faces[ 0 ] = cellID + j;
  faces[ 1 ] = faces[ 0 ] + 1;

  /* The J_DIRECTION faces */
  faces[ 2 ] = cellID + getTotalNumFaces( 0 );
  faces[ 3 ] = faces[ 2 ] + getCellExtent( 0 );

  return 4;
}

//------------------------------------------------------------------------------
inline IndexType StructuredMesh::getCellFacesInternal( IndexType cellID, 
                                                       IndexType j, 
                                                       IndexType k,
                                                       IndexType* faces ) const
{
  SLIC_ASSERT( m_ndims == 3 );
  SLIC_ASSERT( 0 <= cellID && cellID < getNumberOfCells() );
  SLIC_ASSERT( 0 <= j && j < getCellExtent( 1 ) );
  SLIC_ASSERT( 0 <= k && k < getCellExtent( 2 ) );
  SLIC_ASSERT( faces != AXOM_NULLPTR );

  /* The I_DIRECTION faces */
  faces[ 0 ] = cellID + j + getCellExtent( 1 ) * k;
  faces[ 1 ] = faces[ 0 ] + 1;

  /* The J_DIRECTION faces */
  faces[ 2 ] = cellID + getTotalNumFaces( 0 ) + getCellExtent( 0 ) * k;
  faces[ 3 ] = faces[ 2 ] + getCellExtent( 0 );

  /* The K_DIRECTION faces */
  faces[ 4 ] = cellID + getTotalNumFaces( 0 ) + getTotalNumFaces( 1 );
  faces[ 5 ] = faces[ 4 ] + cellKp(); 

  return 6;
}

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_STRUCTUREDMESH_HPP_ */
