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

#ifndef EXTENT_HXX_
#define EXTENT_HXX_

#include "slic/slic.hpp"

// C/C++ includes
#include <cstring> // for memcpy()

namespace axom
{
namespace mint
{

template < typename IndexType >
class Extent
{
public:

  /*!
   * \brief Constructor. Creates an extent instance of the given dimension.
   * \param [in] ndims the number of dimensions.
   * \param [in] ext the extent.
   * \pre ndims >= 1 && ndims <= 3
   */
  Extent( int ndims, const IndexType* ext );

  /*!
   * \brief Returns the dimension of this extent.
   * \return ndims the number of dimensions.
   * \pre ndims >= 1 && ndims <= 3.
   */
  int getDimension() const { return m_ndims; }

  /*!
   * \brief Returns the min along the given dimension.
   * \param [in] idim the dimension in query.
   * \return the min along the given dimension.
   * \pre idim >= 0 && idim < this->getDimension();
   * \pre this->getDimension() >= 1.
   */
  IndexType min( int idim ) const { return m_extent[ idim*2 ]; }

  /*!
   * \brief Returns the max along the given dimension.
   * \param [in] idim the dimension in query.
   * \return the max along the given dimension.
   * \pre idim >= 0 && idim < this->getDimension();
   * \pre this->getDimension() >= 1
   */
  IndexType max( int idim ) const { return m_extent[ idim*2+1 ]; }

  /*!
   * \brief Returns the size along the given dimension, max()-min()+1.
   * \param [in] idim the dimension in query.
   * \return n size of the given dimension.
   */
  int size( int idim ) const
  {
    return( this->max( idim )-this->min( idim )+1 );
  }

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  IndexType jp() const { return m_jp; };

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  IndexType kp() const { return m_kp; };
  /*!
   * \brief Returns the number of nodes covered by this extent instance.
   * \return N the total number of nodes in the extent.
   * \pre this->getDimension() >= 1.
   */
  IndexType getNumNodes() const;

  /*!
   * \brief Returns the number cells covered by this extent instance.
   * \return N the total number of cells in the extent.
   * \pre this->getDimension() >= 1.
   */
  IndexType getNumCells() const;

  /*!
   * \brief Returns the cell offset lookup table
   * \return offsets pointer to the cell offsets table.
   *
   * \note The indent for the cell offsets table
   */
  const IndexType* getCellOffSets() const
  { return &m_cell_offsets[0]; };

  /*!
   * \brief Converts the given grid indices to a one-dimensional linear index.
   * \param [in] i the grid index along the first dimension.
   * \param [in] j the grid index along the second dimension.
   * \param [in] k the grid index along the third dimension.
   * \return linearIdx the linear index.
   * \pre this->getDimension() == 3.
   * \pre i >= 0 && i < this->size( 0 )
   * \pre j >= 0 && j < this->size( 1 )
   * \pre k >= 0 && k < this->size( 2 )
   * \note i,j,k are local grid indices.
   */
  IndexType getLinearIndex( IndexType i, IndexType j, IndexType k ) const;

  /*!
   * \brief Converts the given grid indices to a one-dimensional linear index.
   * \param [in] i the grid index along the first dimension.
   * \param [in] j the grid index along the second dimension.
   * \return linearIdx the linear index.
   * \pre this->getDimension() == 2.
   * \pre i >= 0 && i < this->size( 0 )
   * \pre j >= 0 && j < this->size( 1 )
   * \pre k >= 0 && k < this->size( 2 )
   * \note i,j,k are local grid indices.
   */
  IndexType getLinearIndex( IndexType i, IndexType j ) const;

  /*!
   * \brief Converts the given grid cell indices to a one-dimensional linear
   *  index.
   * \param [in] i the grid cell index along the first dimension.
   * \param [in] j the grid cell index along the second dimension.
   * \param [in] k the grid cell index along the third dimension.
   * \return linearIdx the linear index over the cells.
   * \pre this->getDimension() == 3.
   * \pre i >= 0 && i < this->size( 0 )-1
   * \pre j >= 0 && j < this->size( 1 )-1
   * \pre k >= 0 && k < this->size( 2 )-1
   * \note i,j,k are local grid cell indices.
   */
  IndexType getCellLinearIndex(IndexType i, IndexType j, IndexType k ) const;

  /*!
   * \brief Converts the given grid cell indices to a one-dimensional linear
   *  index.
   * \param [in] i the grid cell index along the first dimension.
   * \param [in] j the grid cell index along the second dimension.
   * \return linearIdx the linear index over the cells.
   * \pre this->getDimension() == 2.
   * \pre i >= 0 && i < this->size( 0 )-1
   * \pre j >= 0 && j < this->size( 1 )-1
   * \note i,j,k are local grid cell indices.
   */
  IndexType getCellLinearIndex( IndexType i, IndexType j ) const;

  /*!
   * \brief Given a one-dimensional linear index, this method computes the
   *  corresponding (i,j) grid indices.
   * \param [in] linearIdx local flat index.
   * \param [out] i the corresponding grid index along the first dimension.
   * \param [out] j the corresponding grid index along the second dimension.
   * \pre this->getDimension() == 2.
   */
  void getGridIndex( IndexType linearIdx, IndexType& i, IndexType& j ) const;

  /*!
   * \brief Given a one-dimensional linear index, this method computes the
   *  corresponding (i,j,k) grid indices.
   * \param [in] linearIdx local flat index.
   * \param [out] i the corresponding grid index along the first dimension.
   * \param [out] j the corresponding grid index along the second dimension.
   * \param [out] k the corresponding grid index along the third dimension.
   * \pre this->getDimension() == 3.
   */
  void getGridIndex( IndexType linearIdx,
                     IndexType& i, IndexType& j, IndexType& k ) const;

private:

  /*!
   * \brief Default constructor.
   * \note Made private to prevent from calling it.
   */
  Extent();

  /*!
   * \brief Builds the cell offsets lookup table.
   * \note Called from the constructor.
   */
  void buildCellOffsets();

  int m_ndims;                  /*!< dimension of this extent    */
  IndexType m_jp;               /*!< stride to the 2nd dimension */
  IndexType m_kp;               /*!< stride to the 3rd dimension */
  IndexType m_extent[6];        /*!< extent of this instance     */
  IndexType m_cell_offsets[8];  /*!< cell offsets                */
};

} /* namespace mint */
} /* namespace axom */

//------------------------------------------------------------------------------
//  Extent Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace mint
{

//------------------------------------------------------------------------------
template < typename IndexType >
Extent< IndexType >::Extent() :
  m_ndims( -1 ),
  m_jp(0),
  m_kp(0)

{
  this->buildCellOffsets();
}

//------------------------------------------------------------------------------
template < typename IndexType >
Extent< IndexType >::Extent( int ndims, const IndexType* ext ) :
  m_ndims( ndims )
{
  SLIC_ASSERT( ndims >= 1 && ndims <= 3 );

  // zero out all extents
  memset(m_extent, 0, 6 * sizeof(IndexType) );

  // copy in the user-supplied extent
  memcpy(m_extent, ext, 2 * ndims * sizeof( IndexType )  );

  // compute strides
  m_jp = 0;
  m_kp = 0;

  if ( ndims > 1 )
  {
    m_jp = this->size( 0 );
  }

  if ( ndims > 2 )
  {
    m_kp = m_jp * this->size( 1 );

  }

  this->buildCellOffsets();
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
void Extent< IndexType >::buildCellOffsets()
{
  m_cell_offsets[ 0 ] = 0;
  m_cell_offsets[ 1 ] = 1;
  m_cell_offsets[ 2 ] = 1 + m_jp;
  m_cell_offsets[ 3 ] = m_jp;

  m_cell_offsets[ 4 ] = m_kp;
  m_cell_offsets[ 5 ] = 1 + m_kp;
  m_cell_offsets[ 6 ] = 1 + m_jp + m_kp;
  m_cell_offsets[ 7 ] = m_jp + m_kp;
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getNumNodes() const
{
  IndexType n = 1;
  for ( int idim=0 ; idim < m_ndims ; ++idim )
  {
    IndexType size = static_cast< IndexType >( this->size( idim ) );
    n *= (size > 0) ? size : 1;
  }
  return( n );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getNumCells() const
{
  IndexType n = 1;
  for ( int idim=0 ; idim < m_ndims ; ++idim )
  {
    IndexType size = static_cast< IndexType >( this->size( idim )-1 );
    n *= ( size > 0) ? size : 1;
  }
  return( n );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getLinearIndex(
  IndexType i, IndexType j, IndexType k ) const
{
  IndexType index = i + j * m_jp + k * m_kp;
  return( index );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getLinearIndex( IndexType i, IndexType j ) const
{
  IndexType index = i + j * m_jp;
  return( index );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getCellLinearIndex(
  IndexType i, IndexType j, IndexType k) const
{
  IndexType cell_jp = (size(0)-1);
  IndexType cell_kp = getDimension() == 3 ? cell_jp * (size(1)-1) : 0;
  IndexType index = i + j * cell_jp  + k * cell_kp;
  return( index );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
IndexType Extent< IndexType >::getCellLinearIndex( IndexType i,
                                                   IndexType j ) const
{
  IndexType cell_jp = (size(0)-1);
  IndexType index = i + j * cell_jp;
  return( index );
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
void Extent< IndexType >::getGridIndex(
  IndexType linearIdx, IndexType& i, IndexType& j) const
{
  SLIC_ASSERT( m_ndims == 2 );

  j = linearIdx / m_jp;
  i = linearIdx - j*m_jp;
}

//------------------------------------------------------------------------------
template < typename IndexType >
inline
void Extent< IndexType >::getGridIndex(
  IndexType linearIdx, IndexType& i, IndexType& j, IndexType& k) const
{
  k = (m_kp > 0) ? (linearIdx / m_kp) : 0;
  j = (m_jp > 0) ? (linearIdx - k*m_kp) / m_jp : 0;
  i = linearIdx - k*m_kp - j*m_jp;
}

} /* namespace mint */
} /* namespae axom */

#endif /* EXTENT_HXX_ */
