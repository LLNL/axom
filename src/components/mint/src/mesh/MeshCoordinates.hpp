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

#ifndef MESHCOORDINATES_HXX_
#define MESHCOORDINATES_HXX_

#include "mint/config.hpp"
#include "mint/Array.hpp"

#include "axom/Macros.hpp"

#define X_COORDINATE 0
#define Y_COORDINATE 1
#define Z_COORDINATE 2

namespace axom
{

#ifdef MINT_USE_SIDRE
namespace sidre
{
class Group;
}
#endif

namespace mint
{

/*!
 * \class MeshCoordinates
 *
 * \brief Provides functionality for storing mesh coordinates.
 */
class MeshCoordinates
{
public:

  /*!
   * \brief TODO:
   * \param dimension
   */
  MeshCoordinates( int dimension );

  /*!
   * \brief Creates a MeshCoordinates instance with the given number of points.
   * \param [in] dimension the dimension of the ambient space.
   * \param [in] capacity the total number of points to allocate space for.
   * \param [in] resize_ratio the ratio to resize the coordinate array by when
   *  it exceeds the capacity. A ratio smaller than 1 prevents dynamic resizing.
   */
  MeshCoordinates( int dimension, IndexType capacity, IndexType size,
                   double resize_ratio );

#ifdef MINT_USE_SIDRE
  /*!
   * \brief Creates a MeshCoordinates instance from a sidre::Group that already
   *  has data.
   * \param [in] group the sidre::Group to use.
   * \pre group != AXOM_NULLPTR.
   */
  MeshCoordinates( sidre::Group* group, int dimension, IndexType size,
                   double resize_ratio );

  /*!
   * \brief Creates a MeshCoordinates instance from an empty sidre::Group.
   * \param [in] group the sidre::Group to use.
   * \param [in] dimension the dimension of the coordinate system.
   * \param [in] capacity the total number of points to allocate space for.
   * \param [in] resize_ratio the ratio to resize the coordinate array by when
   *  it exceeds the capacity. A ratio smaller than 1 prevents dynamic resizing.
   * \pre group != AXOM_NULLPTR.
   */
  MeshCoordinates( sidre::Group* group, int dimension, IndexType capacity,
                   IndexType size, double resize_ratio );
#endif

  /*!
   * \brief Destructor, free's the allocated vectors.
   */
  ~MeshCoordinates();

  /*!
   * \brief Adds a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \pre m_ndims == 1.
   */
  inline void addPoint( double x );

  /*!
   * \brief Adds new points into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] n the number of points to add.
   * \pre m_ndims == 1.
   */
  inline void addPoints( double* x, IndexType n );

  /*!
   * \brief Adds a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \pre m_ndims == 2.
   */
  inline void addPoint( double x, double y );

  /*!
   * \brief Adds new points into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \param [in] n the number of points to add.
   * \pre m_ndims == 2.
   */
  inline void addPoints( double* x, double* y, IndexType n );

  /*!
   * \brief Adds a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \param [in] z the z--coordinate.
   * \pre m_ndims == 3.
   */
  inline void addPoint( double x, double y, double z );

  /*!
   * \brief Adds new points into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \param [in] z the z--coordinate.
   * \param [in] n the number of points to add.
   * \pre m_ndims == 3.
   */
  inline void addPoints( double* x, double* y, double* z, IndexType n );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 1.
   */
  inline void setPoint( IndexType pntIdx, double x );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 2.
   */
  inline void setPoint( IndexType pntIdx, double x, double y );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \param [in] z the z--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 3.
   */
  inline void setPoint( IndexType pntIdx, double x, double y, double z );

  /*!
   * \brief Returns the coordinate of a point at the given dimension.
   * \param [in] pntIdx the index of the point in query.
   * \param [in] dim the dimension in query.
   * \return coord the coordinate of the point.
   * \pre dim < m_ndims
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   */
  inline double getCoordinate( IndexType pntIdx, int dim );

  /*!
   * \brief Returns a const pointer to the coordinate array.
   * \param [in] dim the requested dimension.
   * \return coord_array const pointer to the coordinate array
   * \pre dim < m_ndims
   * \post coord_array != AXOM_NULLPTR.
   */
  inline double* getCoordinateArray( int dim );

  /*!
   * \brief Returns a const pointer to the coordinate array.
   * \param [in] dim the requested dimension.
   * \return coord_array const pointer to the coordinate array
   * \pre dim < m_ndims
   * \post coord_array != AXOM_NULLPTR.
   */
  inline const double* getCoordinateArray( int dim ) const;

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   */
  inline IndexType getCapacity() const
  { return m_coordinates[0]->capacity(); }

  /*!
   * \brief Change the maximum number of points that can currently be held.
   */
  inline void reserve( IndexType capacity );

  /*!
   * \brief Returns the number of points in this MeshCoordinates instance.
   * \return npoint the number points in this MeshCoordinates instance.
   */
  inline IndexType size() const
  { return m_coordinates[0]->size(); }

  /*!
   * \brief Sets the number of points in this MeshCoordinates instance.
   */
  inline void resize( IndexType size );


  inline double getResizeRatio() const
  { return m_coordinates[0]->getResizeRatio(); }


  inline void setResizeRatio( double ratio );

private:
  // TODO: support different memory layouts...

  int m_ndims;
  Array< double >* m_coordinates[ 3 ];

  DISABLE_COPY_AND_ASSIGNMENT(MeshCoordinates);
  DISABLE_MOVE_AND_ASSIGNMENT(MeshCoordinates);

};

//------------------------------------------------------------------------------
// In-lined method implementations
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x )
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ]->append( x );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, IndexType n )
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ]->append( x, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x, double y )
{
  SLIC_ASSERT( m_ndims == 2 );
  m_coordinates[ X_COORDINATE ]->append( x );
  m_coordinates[ Y_COORDINATE ]->append( y );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, double* y, IndexType n )
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ]->append( x, n );
  m_coordinates[ Y_COORDINATE ]->append( y, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x, double y, double z )
{
  SLIC_ASSERT( m_ndims == 3 );
  m_coordinates[ X_COORDINATE ]->append( x );
  m_coordinates[ Y_COORDINATE ]->append( y );
  m_coordinates[ Z_COORDINATE ]->append( z );
}


//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, double* y, double* z,
                                        IndexType n )
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ]->append( x, n );
  m_coordinates[ Y_COORDINATE ]->append( y, n );
  m_coordinates[ Z_COORDINATE ]->append( z, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( IndexType pntIdx, double x )
{
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < size() ) );
  SLIC_ASSERT( m_ndims == 1 );
  (*m_coordinates[ X_COORDINATE ])( pntIdx ) = x;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( IndexType pntIdx, double x, double y )
{
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < size() ) );
  SLIC_ASSERT( m_ndims == 2 );
  (*m_coordinates[ X_COORDINATE ])( pntIdx ) = x;
  (*m_coordinates[ Y_COORDINATE ])( pntIdx ) = y;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( IndexType pntIdx, double x, double y,
                                       double z)
{
  SLIC_ASSERT( m_ndims == 3 );
  (*m_coordinates[ X_COORDINATE ])( pntIdx ) = x;
  (*m_coordinates[ Y_COORDINATE ])( pntIdx ) = y;
  (*m_coordinates[ Z_COORDINATE ])( pntIdx ) = z;
}

//------------------------------------------------------------------------------
inline double MeshCoordinates::getCoordinate( IndexType pntIdx, int dim )
{
  SLIC_ASSERT( dim < m_ndims );
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < size() ) );
  return (*m_coordinates[ dim ])( pntIdx );
}

//------------------------------------------------------------------------------
inline double* MeshCoordinates::getCoordinateArray( int dim )
{
  SLIC_ASSERT( dim < m_ndims );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
inline const double* MeshCoordinates::getCoordinateArray( int dim ) const
{
  SLIC_ASSERT( dim < m_ndims );
  return m_coordinates[ dim ]->getData();
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::reserve( IndexType capacity )
{
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_coordinates[ dim ]->reserve( capacity );
  }
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::resize( IndexType size )
{
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_coordinates[ dim ]->resize( size );
  }
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setResizeRatio( double ratio )
{
  for ( int dim = 0 ; dim < m_ndims ; ++dim )
  {
    m_coordinates[ dim ]->setResizeRatio( ratio );
  }
}

} /* namespace mint */
} /* namespace axom */

#endif /* MESHCOORDINATES_HXX_ */
