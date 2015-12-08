/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file MeshCoordinates.cxx
 *
 * \date Sep 12, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */


#include "quest/MeshCoordinates.hpp"

#include <cassert> // for assert()

namespace meshtk
{

MeshCoordinates::MeshCoordinates(): m_ndims(2)
{
  this->initialize();
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates(int dimension) : m_ndims(dimension)
{
  this->initialize();
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates(int dimension,int npoints) : m_ndims(dimension)
{
  this->initialize(npoints);
}

//------------------------------------------------------------------------------
MeshCoordinates::MeshCoordinates(int dimension, int ndims[3]) :
        m_ndims(dimension)
{
  m_coordinates.resize( m_ndims );
  for ( int i=0; i < m_ndims; ++i ) {
    m_coordinates[ i ].resize( ndims[ i ] );
  }

}

//------------------------------------------------------------------------------
MeshCoordinates::~MeshCoordinates()
{
  for ( int i=0; i < m_ndims; ++i ) {
    m_coordinates[ i ].clear();
  }
  m_coordinates.clear();
}

//------------------------------------------------------------------------------
void MeshCoordinates::insertPoint( double x, double y )
{
  assert( m_ndims == 2 );

  m_coordinates[ X_COORDINATE ].push_back( x );
  m_coordinates[ Y_COORDINATE ].push_back( y );

  assert( m_coordinates[ X_COORDINATE ].size() ==
                     m_coordinates[ Y_COORDINATE ].size() );
}

//------------------------------------------------------------------------------
void MeshCoordinates::insertPoint( double x, double y, double z )
{
  assert( m_ndims == 3 );

  m_coordinates[ X_COORDINATE ].push_back( x );
  m_coordinates[ Y_COORDINATE ].push_back( y );
  m_coordinates[ Z_COORDINATE ].push_back( z );

  assert( m_coordinates[ X_COORDINATE ].size() ==
                     m_coordinates[ Y_COORDINATE ].size() );
  assert( m_coordinates[ X_COORDINATE ].size() ==
                     m_coordinates[ Z_COORDINATE ].size() );
}

//------------------------------------------------------------------------------
void MeshCoordinates::setPoint( int pntIdx, double x, double y )
{
  assert( (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints()) );
  assert( m_ndims == 2 );

  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
}

//------------------------------------------------------------------------------
void MeshCoordinates::setPoint( int pntIdx, double x, double y, double z)
{
  assert( (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints()) );
  assert( m_ndims == 3 );

  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
  m_coordinates[ Z_COORDINATE ][ pntIdx ] = z;
}

//------------------------------------------------------------------------------
double MeshCoordinates::getCoordinate( int pntIdx, int dim )
{
  assert( dim < m_ndims );
  assert( (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints()) );
  return m_coordinates[ dim ][ pntIdx ];
}

//------------------------------------------------------------------------------
double* MeshCoordinates::getCoordinateArray(int dim)
{
  assert( dim < m_ndims );
  return &(m_coordinates[ dim ][ 0 ]);
}

//------------------------------------------------------------------------------
int MeshCoordinates::getNumberOfPoints() const
{
  assert( m_ndims >= 1 );
  return m_coordinates[ X_COORDINATE ].size();
}

//------------------------------------------------------------------------------
void MeshCoordinates::initialize()
{
  assert( m_ndims >= 1 );

  m_coordinates.resize( m_ndims );
  for ( int i=0; i < m_ndims; ++i ) {
    m_coordinates[ i ].reserve( 100 );
  }
}

//------------------------------------------------------------------------------
void MeshCoordinates::initialize( int npoints )
{
  assert( m_ndims >= 1 );

  m_coordinates.resize( m_ndims );
  for ( int i=0; i < m_ndims; ++i ) {
    m_coordinates[ i ].resize( npoints );
  }
}

} /* namespace meshtk */
