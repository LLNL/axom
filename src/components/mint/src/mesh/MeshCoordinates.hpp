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

#include "axom/Macros.hpp"
#include "mint/Vector.hpp"

#define X_COORDINATE 0
#define Y_COORDINATE 1
#define Z_COORDINATE 2

namespace axom
{
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
   * \brief Creates a MeshCoordinates instance with the given number of points.
   * \param [in] dimension the dimension of the ambient space.
   * \param [in] capacity the total number of points to allocate space for.
   * \param [in] resize_ratio the ammount to resize the coordinate array by
   *  when it exceeds the capacity.
   * \pre ndims != AXOM_NULLPTR.
   * \note The points are not initialized but, can be set using setPoint
   */
  MeshCoordinates( int dimension=1, int capacity=100, double resize_ratio=2.0 );

  /*!
   * \brief Destructor.
   */
  ~MeshCoordinates()
  {}

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
  inline void addPoints( double* x, int n );

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
  inline void addPoints( double* x, double* y, int n );

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
  inline void addPoints( double* x, double* y, double* z, int n );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 1.
   */
  inline void setPoint( int pntIdx, double x );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 2.
   */
  inline void setPoint( int pntIdx, double x, double y );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \param [in] z the z--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 3.
   */
  inline void setPoint( int pntIdx, double x, double y, double z );

  /*!
   * \brief Returns the coordinate of a point at the given dimension.
   * \param [in] pntIdx the index of the point in query.
   * \param [in] dim the dimension in query.
   * \return coord the coordinate of the point.
   * \pre dim < m_ndims
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   */
  inline double getCoordinate( int pntIdx, int dim );

  /*!
   * \brief Returns a const pointer to the coordinate array.
   * \param [in] dim the requested dimension.
   * \return coord_array const pointer to the coordinate array
   * \pre dim < m_ndims
   * \post coord_array != AXOM_NULLPTR.
   */
  inline double * getCoordinateArray( int dim );

  /*!
   * \brief Get the maximum number of points that can currently be held.
   * \return N the capacity of m_coordinates.
   */
  inline int getCapacity() const
  { return m_coordinates[0].getCapacity(); }

  /*!
   * \brief Change the maximum number of points that can currently be held.
   */
  inline void setCapacity( int capacity );

  /*!
   * \brief Returns the number of points in this MeshCoordinates instance.
   * \return npoint the number points in this MeshCoordinates instance.
   */
  inline int getSize() const
  { return m_coordinates[0].getSize(); }

  /*!
   * \brief Sets the number of points in this MeshCoordinates instance.
   */
  inline void setSize( int size );

private:
  // TODO: support different memory layouts...

  int m_ndims;
  Vector< double, int > m_coordinates[3]; 

  DISABLE_COPY_AND_ASSIGNMENT(MeshCoordinates);
  DISABLE_MOVE_AND_ASSIGNMENT(MeshCoordinates);
};

//------------------------------------------------------------------------------
// In-lined method implimentation
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x ) 
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ].add( x );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, int n ) 
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ].add( x, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x, double y ) 
{
  SLIC_ASSERT( m_ndims == 2 );
  m_coordinates[ X_COORDINATE ].add( x );
  m_coordinates[ Y_COORDINATE ].add( y );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, double* y, int n ) 
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ].add( x, n );
  m_coordinates[ Y_COORDINATE ].add( y, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoint( double x, double y, double z ) 
{
  SLIC_ASSERT( m_ndims == 3 );
  m_coordinates[ X_COORDINATE ].add( x );
  m_coordinates[ Y_COORDINATE ].add( y );
  m_coordinates[ Z_COORDINATE ].add( z );
}


//------------------------------------------------------------------------------
inline void MeshCoordinates::addPoints( double* x, double* y, double* z, int n ) 
{
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ].add( x, n );
  m_coordinates[ Y_COORDINATE ].add( y, n );
  m_coordinates[ Z_COORDINATE ].add( z, n );
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( int pntIdx, double x ) 
{
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < getSize() ) );
  SLIC_ASSERT( m_ndims == 1 );
  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( int pntIdx, double x, double y ) 
{
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < getSize() ) );
  SLIC_ASSERT( m_ndims == 2 );
  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setPoint( int pntIdx, double x, double y, double z)
{
  SLIC_ASSERT( m_ndims == 3 );
  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
  m_coordinates[ Z_COORDINATE ][ pntIdx ] = z;
}

//------------------------------------------------------------------------------
inline double MeshCoordinates::getCoordinate( int pntIdx, int dim )
{
  SLIC_ASSERT( dim < m_ndims );
  SLIC_ASSERT( ( pntIdx >= 0 ) && ( pntIdx < getSize() ) );
  return m_coordinates[ dim ][ pntIdx ];
}

//------------------------------------------------------------------------------
inline double* MeshCoordinates::getCoordinateArray( int dim )
{
  SLIC_ASSERT( dim < m_ndims );
  return m_coordinates[ dim ].getData();
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setCapacity( int capacity ) 
{
  for ( int dim = 0; dim < m_ndims; ++dim ) {
    m_coordinates[ dim ].setCapacity( capacity );
  }
}

//------------------------------------------------------------------------------
inline void MeshCoordinates::setSize( int size ) 
{
  for ( int dim = 0; dim < m_ndims; ++dim ) {
    m_coordinates[ dim ].setSize( size );
  }
}

} /* namespace mint */
} /* namespace axom */

#endif /* MESHCOORDINATES_HXX_ */
