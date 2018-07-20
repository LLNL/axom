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

#ifndef MINT_EXTENT_HPP_
#define MINT_EXTENT_HPP_

#include "axom/mint/config.hpp"  // for compile-time definitions

#include "axom/slic/interface/slic.hpp"    // for SLIC macros

namespace axom
{
namespace mint
{

static constexpr int I_DIRECTION = 0;
static constexpr int J_DIRECTION = 1;
static constexpr int K_DIRECTION = 2;

/*!
 * \class Extent
 *
 * \brief Defines the logical rectangle corresponding to the regular topology
 *  of a StructuredMesh and provides methods for mapping between the i-j-k
 *  grid index, in the regular topology, to the corresponding linear (flat)
 *  index used to access the arrays that store the field variables and mesh
 *  coordinates.
 *
 * \note This class is used to facilitate operations on structured meshes.
 *
 * \see CurvilinearMesh
 * \see RectilinearMesh
 * \see UniformMesh
 * \see StructuredMesh
 */
class Extent
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  Extent( ) = delete;

  /*!
   * \brief Constructor. Creates an extent instance of the given dimension.
   * \param [in] ndims the number of dimensions.
   * \param [in] ext pointer to buffer consisting of the extent information.
   *
   * \note The supplied `ext` pointer must point to a buffer that has at least
   * \f$ 2 \times N \f$ entries, where N is the dimension of the uniform mesh
   * given in the following order: [imin, imax, jmin, jmax, kmin, kmax]
   *
   * \pre ndims >= 1 && ndims <= 3
   * \pre ext != nullptr
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
  inline const int64& min( int dim ) const
  { return m_extent[ dim*2 ]; }

  /*!
   * \brief Returns the max along the given dimension.
   * \param [in] dim the dimension in query.
   * \return the max along the given dimension.
   * \pre dim >= 0 && dim < getDimension();
   * \pre getDimension() >= 1
   */
  inline const int64& max( int dim ) const
  { return m_extent[ dim*2+1 ]; }

  /*!
   * \brief Returns the size along the given dimension, max()-min()+1.
   * \param [in] dim the dimension in query.
   * \return n size of the given dimension.
   */
  inline IndexType size( int dim ) const
  { return m_sizes[ dim ]; }

  /*!
   * \brief Returns stride to the second dimension.
   * \return jp stride to the second dimension.
   * \post jp >= 0.
   */
  inline IndexType nodeJp() const
  { return m_node_jp; }

  /*!
   * \brief Returns stride to the third dimension.
   * \return kp stride to the third dimension.
   * \post kp >= 0.
   */
  inline IndexType nodeKp() const
  { return m_node_kp; }

  inline IndexType cellJp() const
  { return m_cell_jp; }

  inline IndexType cellKp() const
  { return m_cell_kp; }
  
  /*!
   * \brief Returns the number of nodes covered by this extent instance.
   */
  inline IndexType getNumNodes() const
  { return m_numnodes; }

  /*!
   * \brief Returns the number cells covered by this extent instance.
   */
  inline IndexType getNumCells() const
  { return m_numcells; }

  /*!
   * \brief Returns the number of faces covered by this extent instance.
   */
  inline IndexType getNumFaces() const
  { return m_numfaces; }

  /*!
   * \brief Returns the number edges covered by this extent instance.
   */
  inline IndexType getNumEdges() const
  { return m_numedges; }

  /*!
   * \brief Returns the cell offset lookup table
   * \return offsets pointer to the cell offsets table.
   *
   * \note The indent for the cell offsets table
   */
  inline const IndexType* getCellOffSets() const
  { return &m_cell_offsets[0]; };

  /*!
   * \brief Shifts the given global IJK coordinates of a node or cell to the
   *  local frame of reference of this extent starting from zero, i.e.,
   *  on the interval \f$ [0, N] \f$, where N is the number of nodes (or cells)
   *  along the given dimension.
   *
   * \param [in]  gijk pointer to the global IJK coordinates
   * \param [out] lijk pointer to buffer where the computed local coordinates
   *  will be stored
   *
   * \note The supplied pointers, gijk and lijk must point to buffers that can
   *  hold at least NDIMS values.
   *
   * \pre gijk != nullptr
   * \pre lijk != nullptr
   */
  inline void shiftToLocal( const int64* gijk, IndexType* lijk ) const;

  /*!
   * \brief Shifts the given global IJK coordinates of a node or cell to the
   *  global frame of reference of this extent, i.e., on the interval
   *  [ min( i ), max( i ) ] along each dimension.
   *
   * \param [in]  lijk pointer to buffer that holds the local IJK coordinates
   * \param [out] gijk pointer to buffer where the computed global coordinates
   *  will be stored.
   *
   * \note The supplied pointers, gijk and lijk must point to buffers that can
   *  hold at least NDIMS values.
   *
   * \pre gijk != nullptr
   * \pre lijk != nullptr
   */
  inline void shiftToGlobal( const IndexType* lijk, int64* gijk ) const;

  /*!
   * \brief Converts the given grid indices to a one-dimensional linear index.
   *
   * \param [in] i the grid index along the I_DIRECTION.
   * \param [in] j the grid index along the J_DIRECTION.
   * \param [in] k the grid index along the K_DIRECTION.
   *
   * \return linearIdx the linear index
   *
   * \note The provided i-j-k grid indices are expected to be the shifted
   *  topological coordinates within the local frame of reference for the
   *  extent.
   *
   * \pre getDimension() >= 1.
   * \pre i >= 0 && i < size( 0 )
   * \pre j >= 0 && j < size( 1 ) if getDimesnion() >= 2
   * \pre k >= 0 && k < size( 2 ) if getDimension() == 3
   *
   * \post linearIdx >= 0 && linearIdx < getNumNodes()
   */
  /// @{
  inline IndexType getNodeLinearIndex( const IndexType& i,
                                       const IndexType& j ) const;

  inline IndexType getNodeLinearIndex( const IndexType& i,
                                       const IndexType& j,
                                       const IndexType& k ) const;
  /// @}

  /*!
   * \brief Converts the given cell grid index to a one-dimensional linear
   * index.
   *
   * \param [in] i the grid cell index along the I_DIRECTION.
   * \param [in] j the grid cell index along the J_DIRECTION.
   * \param [in] k the grid cell index along the K_DIRECTION.
   *
   * \return linearIdx the linear index over the cells.
   *
   * \note The provided i-j-k grid indices are expected to be the shifted
   *  topological coordinates within the local frame of reference for the
   *  extent.
   *
   * \pre getDimension() >= 1.
   * \pre i >= 0 && i < size( 0 )-1
   * \pre j >= 0 && j < size( 1 )-1 if getDimension() >= 2
   * \pre k >= 0 && k < size( 2 )-1 if getDimension() == 3
   *
   * \post linearIdx >= 0 && linearIdx < getNumCells()
   */
  /// @{
  inline IndexType getCellLinearIndex( const IndexType& i,
                                       const IndexType& j ) const;

  inline IndexType getCellLinearIndex( const IndexType& i,
                                       const IndexType& j,
                                       const IndexType& k) const;
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

  inline void getNodeGridIndex( const IndexType& linearIdx,
                            IndexType& i,
                            IndexType& j ) const;

  inline void getNodeGridIndex( const IndexType& linearIdx,
                            IndexType& i,
                            IndexType& j,
                            IndexType& k ) const;
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
  inline void getCellGridIndex( const IndexType& linearIdx,
                                IndexType& i,
                                IndexType& j ) const;


  inline void getCellGridIndex( const IndexType& linearIdx,
                                IndexType& i,
                                IndexType& j,
                                IndexType& k ) const;
  /// @}

private:

  /*!
   * \brief Calculate the number of nodes in the mesh given the size of each
   *  dimension.
   */
  void calculateNumberOfNodes();

  /*!
   * \brief Calculate the number of cells in the mesh given the size of each
   *  dimension.
   */
  void calculateNumberOfCells();

  /*!
   * \brief Calculate the number of faces in the mesh given the size of each
   *  dimension.
   */
  void calculateNumberOfFaces();

  /*!
   * \brief Calculate the number of edges in the mesh given the size of each
   *  dimension.
   */
  void calculateNumberOfEdges();

  /*!
   * \brief Builds the cell offsets lookup table.
   * \note Called from the constructor.
   */
  void buildCellOffsets();

  int m_ndims;                                // dimension of this extent   
  IndexType m_numnodes;                       // the number of nodes        
  IndexType m_numcells;                       // the number of cells        
  IndexType m_numfaces;                       // the number of faces        
  IndexType m_numedges;                       // the number of edges        
  IndexType m_node_jp;                        // stride to the 2nd dimension
  IndexType m_node_kp;                        // stride to the 3rd dimension
  IndexType m_cell_jp;                        //
  IndexType m_cell_kp;                        //
  IndexType m_sizes[3] = { 1, 1, 1 };         // size along each dimension  
  int64 m_extent[6] = { 0, 0, 0, 0, 0, 0 };   // extent of this instance    
  IndexType m_cell_offsets[ 8 ];              // cell offsets               
};


//------------------------------------------------------------------------------
//  In-lined Extent Implementation
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void Extent::shiftToLocal( const int64* gijk, IndexType* lijk ) const
{
  SLIC_ASSERT( gijk != nullptr );
  SLIC_ASSERT( lijk != nullptr );

  const int N = getDimension( );
  for ( int i=0 ; i < N ; ++i )
  {
    lijk[ i ] = static_cast< IndexType >( gijk[ i ] - min( i ) );
  }

}

//------------------------------------------------------------------------------
inline void Extent::shiftToGlobal( const IndexType* lijk, int64* gijk ) const
{
  SLIC_ASSERT( lijk != nullptr );
  SLIC_ASSERT( gijk != nullptr );

  const int N = getDimension();
  for ( int i=0 ; i < N ; ++i )
  {
    gijk[ i ] = lijk[ i ] + min( i );
  }
}

//------------------------------------------------------------------------------
inline IndexType Extent::getNodeLinearIndex( const IndexType& i,
                                             const IndexType& j,
                                             const IndexType& k ) const
{
  SLIC_ASSERT( m_ndims==3 );
  return i + j*m_node_jp + k*m_node_kp;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getNodeLinearIndex( const IndexType& i,
                                             const IndexType& j ) const
{
  SLIC_ASSERT( m_ndims==2 );
  return i + j*m_node_jp;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getCellLinearIndex( const IndexType& i,
                                             const IndexType& j,
                                             const IndexType& k ) const
{
  SLIC_ASSERT( m_ndims==3 );
  return i + j * m_cell_jp  + k * m_cell_kp;
}

//------------------------------------------------------------------------------
inline IndexType Extent::getCellLinearIndex( const IndexType& i,
                                             const IndexType& j ) const
{
  SLIC_ASSERT( m_ndims==2 );
  return i + j * m_cell_jp;
}

//------------------------------------------------------------------------------
inline void Extent::getNodeGridIndex( const IndexType& linearIdx,
                                      IndexType& i,
                                      IndexType& j) const
{
  SLIC_ASSERT( m_ndims==2 );
  j = linearIdx / m_node_jp;
  i = linearIdx - j*m_node_jp;
}

//------------------------------------------------------------------------------
inline void Extent::getNodeGridIndex( const IndexType& linearIdx,
                                      IndexType& i,
                                      IndexType& j,
                                      IndexType& k) const
{
  SLIC_ASSERT( m_ndims==3 );
  k = linearIdx / m_node_kp;
  j = (linearIdx - k*m_node_kp) / m_node_jp;
  i = linearIdx - k*m_node_kp - j*m_node_jp;
}

//------------------------------------------------------------------------------
inline void Extent::getCellGridIndex( const IndexType& linearIdx,
                                      IndexType& i,
                                      IndexType& j ) const
{
  SLIC_ASSERT( m_ndims==2 );
  i = linearIdx % m_cell_jp;
  j = linearIdx / m_cell_jp;
}

//------------------------------------------------------------------------------
inline void Extent::getCellGridIndex( const IndexType& linearIdx,
                                      IndexType& i,
                                      IndexType& j,
                                      IndexType& k ) const
{
  SLIC_ASSERT( m_ndims==3 );
  k = linearIdx / m_cell_kp;
  j = (linearIdx - k*m_cell_kp) / m_cell_jp;
  i = linearIdx - k*m_cell_kp - j*m_cell_jp;
}


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_EXTENT_HPP_ */
