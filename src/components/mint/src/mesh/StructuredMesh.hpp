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

#include "axom/Types.hpp"       // for axom types
#include "axom/Macros.hpp"      // for axom macros

#include "mint/config.hpp"      // for compile-time definitions
#include "mint/CellTypes.hpp"   // for the CellTypes enum and definitions
#include "mint/Extent.hpp"      // for mint::Extent
#include "mint/Mesh.hpp"        // for mint::Mesh base class

#include "slic/slic.hpp"        // for SLIC macros

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
  virtual ~StructuredMesh();

/// \name Cells
/// @{

  /*!
   * \brief Return the number of cells in the mesh.
   */
  virtual IndexType getNumberOfCells() const final override
  { return m_extent->getNumCells(); }

  /*!
   * \brief Return the number of nodes associated with the given cell.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   *
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType
  getNumberOfCellNodes(IndexType AXOM_NOT_USED(cellID)=0 ) const override final;

  /*!
   * \brief Return the type of cell this mesh holds. SEGMENT, QUAD, or HEX
   *  depending on the dimension.
   *
   * \param [in] cellID the ID of the cell in question, this parameter is
   *  ignored.
   */
  /// @{

  virtual CellType getCellType( IndexType cellID=0 ) const override final;

  /// @}

  /*!
   * \brief Copy the connectivity of the given cell into the provided buffer.
   *  The buffer must be of length at least getNumberOfCellNodes( cellID ).
   *
   * \param [in] cellID the ID of the cell in question.
   * \param [out] cell the buffer into which the connectivity is copied, must
   *  be of length at least getNumberOfCellNodes().
   *
   * \return The number of nodes for the given cell.
   *
   * \pre cell != AXOM_NULLPTR
   * \pre 0 <= cellID < getNumberOfCells()
   */
  virtual IndexType
  getCell( IndexType cellID, IndexType* cell ) const override final;

/// @}

/// \name Nodes
/// @{

  /*!
   * \brief Return the number of nodes in the mesh.
   */
  virtual IndexType getNumberOfNodes() const final override
  { return m_extent->getNumNodes(); }

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
  virtual void getNode( IndexType nodeID, double* node ) const override = 0;

  /*!
   * \brief Returns pointer to the requested mesh coordinate buffer.
   *
   * \param [in] dim the dimension of the requested coordinate buffer
   * \return ptr pointer to the coordinate buffer.
   *
   * \note if hasExplicitCoordinates() == true then the length of the returned
   *  buffer is getNumberOfNodes(). Otherwise the UniformMesh returns
   *  AXOM_NULLPTR and the RectilinearMesh returns a pointer to the associated
   *  dimension scale which is of length
   *  static_cast< RectilinearMesh* >( this )->getNumberOfNodesAlongDim( dim ).
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
  {
    SLIC_ERROR( "NOT IMPLEMENTED!!!" );
    return 0;
  }

/// @}

/// \name Edges
/// @{

  /*!
   * \brief Return the number of edges in the mesh.
   */
  virtual IndexType getNumberOfEdges() const final override
  {
    SLIC_ERROR( "NOT IMPLEMENTED!!!" );
    return 0;
  }

/// @}

  virtual bool isExternal() const final override
  { return false; }

/// @}

/// \name Attribute Querying Methods
/// @{

  /*!
   * \brief Returns the dimensions of this mesh instance.
   * \param [in] ndims 3-tuple holding the number of nodes in each dimension.
   * \post ndims[ i ] >= 1, \$ \forall i \in [0,2] \$
   */
  inline void getExtentSize( IndexType ndims[ 3 ] ) const
  {
    for ( int i=0 ; i < 3 ; ++i )
    {
      ndims[ i ] = m_extent->size( i );
    } // END for
  }

  /*!
   * \brief Returns the number of nodes along the given dimension.
   * \param [in] dim the dimension to querry.
   * \pre 0 <= dim < 3
   */
  inline IndexType getNumberOfNodesAlongDim( IndexType dim ) const
  {
    SLIC_ASSERT( 0 <= dim && dim < 3 );
    return m_extent->size( dim );
  }

/// @}

/// \name Indexing Helper Methods
/// @{

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  inline IndexType jp() const
  { return m_extent->jp(); }

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  inline IndexType kp() const
  { return m_extent->kp(); }

  /*!
   * \brief Gets a const reference to the Extent object of this StructuredMesh
   * \return ext const reference to the Extent object
   *
   * \post ext != AXOM_NULLPTR
   * \see Extent
   */
  inline const Extent* getExtent() const
  { return m_extent; };

  /*!
   * \brief Returns the linear index corresponding to the given logical node
   *  indices.
   * \param [in] i logical node index of the first dimension.
   * \param [in] j logical node index of the second dimension.
   * \param [in] k logical node index of the third dimension (optional)
   * \return idx the corresponding linear index of the node.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   */
  /// @{

  inline IndexType
  getLinearIndex( IndexType i, IndexType j, IndexType k  ) const
  { return m_extent->getLinearIndex( i, j, k ); };

  inline IndexType getLinearIndex( IndexType i, IndexType j ) const
  { return m_extent->getLinearIndex( i, j ); };

  /// @}

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   * \return idx the corresponding linear index of the cell.
   *
   * \note Each method is valid only for the appropriate dimension of the mesh.
   */
  /// @{

  inline IndexType
  getCellLinearIndex( IndexType i, IndexType j, IndexType k ) const
  { return m_extent->getCellLinearIndex( i, j, k); };

  inline IndexType getCellLinearIndex( IndexType i, IndexType j ) const
  { return m_extent->getCellLinearIndex( i, j ); };

  /// @}

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
  inline void getCell( IndexType i, IndexType j, IndexType* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j,k)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre getDimension() == 3.
   */
  inline void
  getCell( IndexType i, IndexType j, IndexType k, IndexType* cell) const;

/// @}

/// @}

protected:

  /*!
   * \brief Constructs a structured mesh instance from the given extent.
   * \param [in] meshType the mesh type
   * \param [in] dimension the mesh dimension
   * \param [in] ext the structured mesh extent.
   */
  StructuredMesh( int meshType, int dimension, const int64* ext );

  /*!
   * \brief Creates a structured mesh of specified type and dimension
   * \param [in] meshType the mesh type
   * \param [in] dimension the mesh dimension
   */
  StructuredMesh( int meshType, int dimension );

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
   * \pre group !=  AXOM_NULLPTR
   * \pre blueprint::validRootGroup( group ) == true
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
   * \pre group != AXOM_NULLPTR
   * \pre group->getNumGroups() == 0
   * \pre group->getNumViews() == 0
   *
   * \post blueprint::validRootGroup( group )
   *
   * \note This constructor forwards this call to the parent Mesh class.
   *
   * \see Mesh( int ndims, int type, sidre::Group*,
   *            const std::string& topo, const std::string& coordset );
   */
  StructuredMesh( int meshType, int dimension,
                  sidre::Group* group,
                  const std::string& topo,
                  const std::string& coordset );


#endif

  /*!
   * \brief Initializes members of the mesh class.
   */
  void init( );

  /*!
   * \brief Helper method to allocate FieldData on the mesh.
   */
  void allocateFields( );

  Extent* m_extent; /*!< grid extent */

private:
  DISABLE_COPY_AND_ASSIGNMENT( StructuredMesh );
  DISABLE_MOVE_AND_ASSIGNMENT( StructuredMesh );
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline int
StructuredMesh::getNumberOfCellNodes( IndexType AXOM_NOT_USED(cellID) ) const
{
  const int cell_type = getCellType();
  return cell_info[ cell_type ].num_nodes;
}

//------------------------------------------------------------------------------
inline
IndexType StructuredMesh::getCell( IndexType cellID, IndexType* cell ) const
{
  SLIC_ASSERT( cell != AXOM_NULLPTR );
  SLIC_ASSERT( 0 <= cellID && cellID < getNumberOfCells() );

  const IndexType* offsets_table = m_extent->getCellOffSets();
  const IndexType num_cell_nodes = getNumberOfCellNodes();

  // Calculate logical indices of the cell's first corner node.
  IndexType ii, jj, kk;
  m_extent->getCellGridIndex( cellID, ii, jj, kk );

  // Use the offsets table to get the all the cell nodes.
  const IndexType n0 = getLinearIndex( ii, jj, kk );
  for ( IndexType i = 0 ; i < num_cell_nodes; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  }
  return num_cell_nodes;
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell( IndexType i, IndexType j,
                                     IndexType* cell ) const
{
  SLIC_ASSERT( getDimension() == 2 );
  SLIC_ASSERT( cell != AXOM_NULLPTR );

  const IndexType* offsets_table = m_extent->getCellOffSets();
  const IndexType n0 = m_extent->getLinearIndex(i, j);

  for ( IndexType i = 0 ; i < 4 ; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  }

}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell( IndexType i, IndexType j, IndexType k,
                                     IndexType* cell ) const
{
  SLIC_ASSERT( getDimension() == 3 );
  SLIC_ASSERT( cell != AXOM_NULLPTR );

  const IndexType* offsets_table = m_extent->getCellOffSets();
  const IndexType n0 = m_extent->getLinearIndex(i, j, k);

  for ( IndexType i = 0 ; i < 8 ; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  }
}

//------------------------------------------------------------------------------
inline
CellType StructuredMesh::getCellType( IndexType AXOM_NOT_USED(cellID) ) const
{
  return ( (m_ndims==3)? mint::HEX :
             ( (m_ndims==2)? mint::QUAD : mint::SEGMENT ) );
}

} /* namespace mint */
} /* namespace axom */

#endif /* STRUCTUREDMESH_HXX_ */
