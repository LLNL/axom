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

#ifndef STRUCTUREDMESH_HXX_
#define STRUCTUREDMESH_HXX_

#include "axom/Types.hpp"
#include "axom/Macros.hpp"

#include "mint/CellTypes.hpp"
#include "mint/config.hpp"
#include "mint/Extent.hpp"
#include "mint/Mesh.hpp"

#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

class StructuredMesh : public Mesh
{
public:

  /*!
   * \brief Destructor.
   */
  virtual ~StructuredMesh() { }

  /// \name Structured Mesh API
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
      ndims[ i ] = m_extent.size( i );
    } // END for

  }

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  inline IndexType jp() const
  { return m_extent.jp(); }

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  inline IndexType kp() const
  { return m_extent.kp(); }


  /*!
   * \brief Returns the number of cell nodes.
   * \return N the number of nodes per cell.
   * \post N=4 if 2-D, else N=8.
   */
  inline int getNumberOfCellNodes() const;

  /*!
   * \brief Returns the cell type.
   * \return ctype the cell type
   * \see CellType.hpp
   */
  inline int getCellType( ) const;

  /*!
   * \brief Returns the linear index corresponding to the given logical indices.
   * \param [in] i logical index of the first dimension.
   * \param [in] j logical index of the second dimension.
   * \param [in] k logical index of the third dimension (optional)
   * \return idx the corresponding linear index.
   */
  inline IndexType getLinearIndex( IndexType i, IndexType j,
                                    IndexType k  ) const
  { return m_extent.getLinearIndex( i, j, k ); };

  /*!
   * \brief Returns the linear index corresponding to the given logical indices.
   * \param [in] i logical index of the first dimension.
   * \param [in] j logical index of the second dimension.
   * \return idx the corresponding linear index.
   */
  inline IndexType getLinearIndex( IndexType i, IndexType j ) const
  { return m_extent.getLinearIndex( i, j ); };

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell
   * indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   * \return idx the corresponding linear index of the cell.
   */
  inline IndexType getCellLinearIndex( IndexType i, IndexType j,
                                        IndexType k ) const
  { return m_extent.getCellLinearIndex( i, j, k); };

  /*!
   * \brief Returns the linear index corresponding to the given logical grid
   * cell
   * indices.
   * \param [in] i logical cell index of the first dimension.
   * \param [in] j logical cell index of the second dimension.
   * \param [in] k logical cell index of the third dimension (optional)
   * \return idx the corresponding linear index of the cell.
   */
  inline IndexType getCellLinearIndex( IndexType i, IndexType j ) const
  { return m_extent.getCellLinearIndex( i, j ); };

  /*!
   * \brief Returns the cell connectivity for the given cell.
   * \param [in] cellIdx the index of the cell in query.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre cellIdx >= 0 && cellIdx < getNumberOfCells()
   * \pre the user-supplied cell buffer must be of getNumberOfCellNodes() size.
   */
  inline void getCell( IndexType cellIdx, IndexType* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre this->getDimension() == 2.
   */
  inline void getCell( IndexType i, IndexType j, IndexType* cell ) const;

  /*!
   * \brief Returns the cell connectivity of the cell at (i,j,k)
   * \param [in] i logical index of the cell along the first dimension.
   * \param [in] j logical index of the cell along the second dimension.
   * \param [in] k logical index of the cell along the third dimension.
   * \param [out] cell pointer to buffer to populate with the cell connectivity.
   * \pre this->getDimension() == 3.
   */
  inline void getCell( IndexType i, IndexType j, IndexType k,
                       IndexType* cell) const;


  /// @}

protected:

  /*!
   * \brief Default constructor.
   * \note Made private to prevent users from calling it.
   */
  StructuredMesh();

  /*!
   * \brief Constructs a structured mesh instance from the given extent.
   * \param [in] ext the structured mesh extent.
   */
  StructuredMesh( int meshType, int ndims, const int64 ext[6] );

  /*!
   * \brief Constructs a structured mesh instance from the given extent that is
   *  identified by the given blockId and partitionId pair.
   * \param [in] meshType the structured mesh type.
   * \param [in] ext the structured mesh extent.
   * \param [in] blockId the block ID of the mesh.
   * \param [in] partId the partition ID of the mesh.
   */
  StructuredMesh( int meshType, int ndims, const int64 ext[6],
                  int blockId,
                  int partId );

  Extent m_extent; /*!< grid extent */

private:
  DISABLE_COPY_AND_ASSIGNMENT( StructuredMesh );
  DISABLE_MOVE_AND_ASSIGNMENT( StructuredMesh );
};


//------------------------------------------------------------------------------
//      In-lined Method Implementations
//------------------------------------------------------------------------------

inline int StructuredMesh::getCellType( ) const
{
  return ( (m_ndims==3)? mint::HEX :
             ( (m_ndims==2)? mint::QUAD : mint::SEGMENT ) );
}

//------------------------------------------------------------------------------
inline int StructuredMesh::getNumberOfCellNodes( ) const
{
  const int cell_type = getCellType( );
  return cell_info[ cell_type ].num_nodes;
}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell( IndexType cellIdx,
                                     IndexType* cell ) const
{
  SLIC_ASSERT( cell != AXOM_NULLPTR );
  SLIC_ASSERT( (cellIdx >= 0) && (cellIdx < getNumberOfCells() ) );

  const IndexType* offsets_table = m_extent.getCellOffSets();
  const IndexType num_cell_nodes = getNumberOfCellNodes();

  // STEP 0: calculate logical indices of the cell's first corner node.
  const IndexType jp_minus_1 = jp() - 1;
  IndexType ii = 0, jj = 0, kk = 0;

  if  ( this->getDimension() == 1 )
  {
    SLIC_ASSERT( num_cell_nodes == 2 );
    ii = cellIdx;
  }
  else if ( this->getDimension() == 2 )
  {
    SLIC_ASSERT( num_cell_nodes == 4 );

    ii = cellIdx % jp_minus_1;
    jj = cellIdx / jp_minus_1;
  }
  else
  {
    SLIC_ASSERT(  this->getDimension() == 3 );
    SLIC_ASSERT(  num_cell_nodes == 8 );

    IndexType kp_minus_1 = m_extent.size(1)-1;
    ii = cellIdx % jp_minus_1;
    jj = cellIdx / jp_minus_1 % kp_minus_1;
    kk = cellIdx / ( jp_minus_1 * kp_minus_1 );
  }

  // STEP 1: calculate linear index of corner node from (ii,jj,kk)
  const IndexType n0  = ii + jj*m_extent.jp() + kk*m_extent.kp();

  // STEP 2: Last, use the offsets table to get the all the cell nodes.
  for ( IndexType i=0 ; i < num_cell_nodes ; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  } // END for all cell nodes

}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell( IndexType i, IndexType j,
                                     IndexType* cell ) const
{
  SLIC_ASSERT( this->getDimension() == 2 );
  SLIC_ASSERT( cell != AXOM_NULLPTR );

  const IndexType* offsets_table = m_extent.getCellOffSets();
  const IndexType n0 = m_extent.getLinearIndex(i, j);

  for ( IndexType i = 0 ; i < 4 ; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  }

}

//------------------------------------------------------------------------------
inline void StructuredMesh::getCell( IndexType i, IndexType j, IndexType k,
                                     IndexType* cell ) const
{
  SLIC_ASSERT( this->getDimension() == 3 );
  SLIC_ASSERT( cell != AXOM_NULLPTR );

  const IndexType* offsets_table = m_extent.getCellOffSets();
  const IndexType n0 = m_extent.getLinearIndex(i, j, k);

  for ( IndexType i = 0 ; i < 8 ; ++i )
  {
    cell[ i ] = n0 + offsets_table[ i ];
  }
}


} /* namespace mint */
} /* namespace axom */

#endif /* STRUCTUREDMESH_HXX_ */
