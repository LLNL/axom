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
#include "mint/config.hpp"

// C/C++ includes
#include <cstring> // for memcpy()

namespace axom
{
namespace mint
{

class Extent
{
public:

  /*!
   * \brief Constructor. Creates an extent instance of the given dimension.
   * \param [in] ndims the number of dimensions.
   * \param [in] ext the extent.
   * \pre ndims >= 1 && ndims <= 3
   */
  Extent( int ndims, const int64* ext );

  /*!
   * \brief Returns the dimension of this extent.
   * \return ndims the number of dimensions.
   * \pre ndims >= 1 && ndims <= 3.
   */
  inline int getDimension() const
  { return m_ndims; }

  /*!
   * \brief Returns the min along the given dimension.
   * \param [in] dim the dimension in query.
   * \return the min along the given dimension.
   * \pre dim >= 0 && dim < getDimension();
   * \pre getDimension() >= 1.
   */
  inline int64 min( int dim ) const
  { return m_extent[ dim*2 ]; }

  /*!
   * \brief Returns the max along the given dimension.
   * \param [in] dim the dimension in query.
   * \return the max along the given dimension.
   * \pre dim >= 0 && dim < getDimension();
   * \pre getDimension() >= 1
   */
  inline int64 max( int dim ) const
  { return m_extent[ dim*2+1 ]; }

  /*!
   * \brief Returns the size along the given dimension, max()-min()+1.
   * \param [in] dim the dimension in query.
   * \return n size of the given dimension.
   */
  inline IndexType size( int dim ) const
  { return max( dim ) - min( dim ) + 1; }

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  inline IndexType jp() const
  { return m_jp; };

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  inline IndexType kp() const
  { return m_kp; };
  /*!
   * \brief Returns the number of nodes covered by this extent instance.
   * \return N the total number of nodes in the extent.
   * \pre getDimension() >= 1.
   */
  inline IndexType getNumNodes() const;

  /*!
   * \brief Returns the number cells covered by this extent instance.
   * \return N the total number of cells in the extent.
   * \pre getDimension() >= 1.
   */
  inline IndexType getNumCells() const;

  /*!
   * \brief Returns the cell offset lookup table
   * \return offsets pointer to the cell offsets table.
   *
   * \note The indent for the cell offsets table
   */
  inline const IndexType* getCellOffSets() const
  { return &m_cell_offsets[0]; };

  /*!
   * \brief Converts the given grid indices to a one-dimensional linear index.
   * \param [in] i the grid index along the first dimension.
   * \param [in] j the grid index along the second dimension.
   * \param [in] k the grid index along the third dimension.
   * \return linearIdx the linear index.
   * \pre getDimension() == 3.
   * \pre i >= 0 && i < size( 0 )
   * \pre j >= 0 && j < size( 1 )
   * \pre k >= 0 && k < size( 2 )
   * \note i,j,k are local grid indices.
   */
  inline IndexType getLinearIndex( IndexType i, IndexType j,
                                    IndexType k ) const;

  /*!
   * \brief Converts the given grid indices to a one-dimensional linear index.
   * \param [in] i the grid index along the first dimension.
   * \param [in] j the grid index along the second dimension.
   * \return linearIdx the linear index.
   * \pre getDimension() == 2.
   * \pre i >= 0 && i < size( 0 )
   * \pre j >= 0 && j < size( 1 )
   * \pre k >= 0 && k < size( 2 )
   * \note i,j,k are local grid indices.
   */
  inline IndexType getLinearIndex( IndexType i, IndexType j ) const;

  /*!
   * \brief Converts the given grid cell indices to a one-dimensional linear
   * index.
   * \param [in] i the grid cell index along the first dimension.
   * \param [in] j the grid cell index along the second dimension.
   * \param [in] k the grid cell index along the third dimension.
   * \return linearIdx the linear index over the cells.
   * \pre getDimension() == 3.
   * \pre i >= 0 && i < size( 0 )-1
   * \pre j >= 0 && j < size( 1 )-1
   * \pre k >= 0 && k < size( 2 )-1
   * \note i,j,k are local grid cell indices.
   */
  inline IndexType getCellLinearIndex( IndexType i, IndexType j,
                                        IndexType k ) const;

  /*!
   * \brief Converts the given grid cell indices to a one-dimensional linear
   * index.
   * \param [in] i the grid cell index along the first dimension.
   * \param [in] j the grid cell index along the second dimension.
   * \return linearIdx the linear index over the cells.
   * \pre getDimension() == 2.
   * \pre i >= 0 && i < size( 0 )-1
   * \pre j >= 0 && j < size( 1 )-1
   * \note i,j,k are local grid cell indices.
   */
  inline IndexType getCellLinearIndex( IndexType i, IndexType j ) const;

  /*!
   * \brief Given a one-dimensional linear index, this method computes the
   *  corresponding (i,j) grid indices.
   * \param [in] linearIdx local flat index.
   * \param [out] i the corresponding grid index along the first dimension.
   * \param [out] j the corresponding grid index along the second dimension.
   * \pre getDimension() == 2.
   */
  inline void getGridIndex( IndexType linearIdx, IndexType& i,
                            IndexType& j ) const;

  /*!
   * \brief Given a one-dimensional linear index, this method computes the
   *  corresponding (i,j,k) grid indices.
   * \param [in] linearIdx local flat index.
   * \param [out] i the corresponding grid index along the first dimension.
   * \param [out] j the corresponding grid index along the second dimension.
   * \param [out] k the corresponding grid index along the third dimension.
   * \pre getDimension() == 3.
   */
  inline void getGridIndex( IndexType linearIdx,
                            IndexType& i, IndexType& j, IndexType& k ) const;


  inline void getCellGridIndex( IndexType linearIdx, IndexType& i,
                                IndexType& j ) const;

  inline void getCellGridIndex( IndexType linearIdx, IndexType& i, 
                                IndexType& j, IndexType& k ) const;

private:

  /*!
   * \brief Builds the cell offsets lookup table.
   * \note Called from the constructor.
   */
  inline void buildCellOffsets();

  int m_ndims;                   /* dimension of this extent    */
  IndexType m_jp;                /* stride to the 2nd dimension */
  IndexType m_kp;                /* stride to the 3rd dimension */
  int64 m_extent[6];             /* extent of this instance     */
  IndexType m_cell_offsets[8];   /* cell offsets                */

  /*!
   * \brief Default constructor.
   * \note Made private to prevent from calling it.
   */
  Extent();
};


//------------------------------------------------------------------------------
//  In-lined Extent Implementation
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void Extent::buildCellOffsets()
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
inline IndexType Extent::getNumNodes() const
{
  IndexType n = 1;
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    IndexType n_dim = static_cast< IndexType >( size( dim ) );
    n *= (n_dim > 0) ? n_dim : 1;
  }
  return n;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getNumCells() const
{
  IndexType n = 1;
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    IndexType n_dim = static_cast< IndexType >( size( dim ) - 1 );
    n *= ( n_dim > 0) ? n_dim : 1;
  }
  return n;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getLinearIndex( IndexType i, IndexType j,
                                          IndexType k ) const
{ return i + j * m_jp + k * m_kp; }

//------------------------------------------------------------------------------
inline IndexType Extent::getLinearIndex( IndexType i, IndexType j ) const
{ return i + j * m_jp; }

//------------------------------------------------------------------------------
inline IndexType Extent::getCellLinearIndex( IndexType i, IndexType j,
                                              IndexType k) const
{
  const IndexType jp_minus_1 = jp() - 1;
  const IndexType kp_minus_1 = kp() - 1;
  const IndexType cell_jp = jp_minus_1;
  const IndexType cell_kp = getDimension() == 3 ? cell_jp * kp_minus_1 : 0;
  return i + j * cell_jp  + k * cell_kp;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getCellLinearIndex( IndexType i, IndexType j ) const
{
  IndexType cell_jp = (size(0)-1);
  IndexType index = i + j * cell_jp;
  return index;
}

//------------------------------------------------------------------------------
inline void Extent::getGridIndex( IndexType linearIdx, IndexType& i,
                                  IndexType& j) const
{
  SLIC_ASSERT( m_ndims == 2 );

  j = linearIdx / m_jp;
  i = linearIdx - j*m_jp;
}

//------------------------------------------------------------------------------
inline void Extent::getGridIndex( IndexType linearIdx, IndexType& i,
                                  IndexType& j, IndexType& k) const
{
  k = (kp() > 0) ? (linearIdx / kp()) : 0;
  j = (jp() > 0) ? (linearIdx - k*kp()) / jp() : 0;
  i = linearIdx - k*kp() - j*jp();
}

//------------------------------------------------------------------------------
inline void Extent::getCellGridIndex( IndexType linearIdx, IndexType& i,
                                      IndexType& j ) const
{
  SLIC_ASSERT( m_ndims == 2 );
  const IndexType jp_minus_1 = jp() - 1;
  i = linearIdx % jp_minus_1;
  j = linearIdx / jp_minus_1;
}

//------------------------------------------------------------------------------
inline void Extent::getCellGridIndex( IndexType linearIdx, IndexType& i, 
                                      IndexType& j, IndexType& k ) const
{ 
  const IndexType jp_minus_1 = jp() - 1;
  const IndexType kp_minus_1 = kp() - 1;
  k = (kp() > 0) ? linearIdx / ( jp_minus_1 * kp_minus_1 ) : 0;
  linearIdx -= k * kp_minus_1;
  j = (jp() > 0) ? linearIdx / jp_minus_1 : 0;
  linearIdx -= j * jp_minus_1;
  i = linearIdx;
}


} /* namespace mint */
} /* namespae axom */

#endif /* EXTENT_HXX_ */
