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

#include "mint/MeshCoordinates.hpp"

#include "slic/slic.hpp"

namespace axom
{
namespace mint
{

MeshCoordinates::MeshCoordinates() : m_ndims(2)
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
  for ( int i=0 ; i < m_ndims ; ++i )
  {
    m_coordinates[ i ].resize( ndims[ i ] );
  }

}

//------------------------------------------------------------------------------
MeshCoordinates::~MeshCoordinates()
{
  for ( int i=0 ; i < m_ndims ; ++i )
  {
    m_coordinates[ i ].clear();
  }
  m_coordinates.clear();
}

//------------------------------------------------------------------------------
void MeshCoordinates::insertPoint( double x )
{
  SLIC_ASSERT( m_ndims == 1 );

  m_coordinates[ X_COORDINATE ].push_back( x );
}

//------------------------------------------------------------------------------
void MeshCoordinates::insertPoint( double x, double y )
{
  SLIC_ASSERT( m_ndims == 2 );

  m_coordinates[ X_COORDINATE ].push_back( x );
  m_coordinates[ Y_COORDINATE ].push_back( y );

  SLIC_ASSERT( m_coordinates[ X_COORDINATE ].size() ==
               m_coordinates[ Y_COORDINATE ].size() );
}

//------------------------------------------------------------------------------
void MeshCoordinates::insertPoint( double x, double y, double z )
{
  SLIC_ASSERT( m_ndims == 3 );

  m_coordinates[ X_COORDINATE ].push_back( x );
  m_coordinates[ Y_COORDINATE ].push_back( y );
  m_coordinates[ Z_COORDINATE ].push_back( z );

  SLIC_ASSERT(  m_coordinates[ X_COORDINATE ].size() ==
                m_coordinates[ Y_COORDINATE ].size() );
  SLIC_ASSERT(  m_coordinates[ X_COORDINATE ].size() ==
                m_coordinates[ Z_COORDINATE ].size() );
}

//------------------------------------------------------------------------------
void MeshCoordinates::setPoint( int pntIdx, double x )
{
  SLIC_ASSERT(  ( pntIdx >= 0 ) && ( pntIdx < this->getNumberOfPoints() ) );
  SLIC_ASSERT(  m_ndims == 1 );

  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
}

//------------------------------------------------------------------------------
void MeshCoordinates::setPoint( int pntIdx, double x, double y )
{
  SLIC_ASSERT(  ( pntIdx >= 0 ) && ( pntIdx < this->getNumberOfPoints() ) );
  SLIC_ASSERT(  m_ndims == 2 );

  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
}

//------------------------------------------------------------------------------
void MeshCoordinates::setPoint( int pntIdx, double x, double y, double z)
{
  SLIC_ASSERT(  ( pntIdx >= 0 ) && ( pntIdx < this->getNumberOfPoints() ) );
  SLIC_ASSERT(  m_ndims == 3 );

  m_coordinates[ X_COORDINATE ][ pntIdx ] = x;
  m_coordinates[ Y_COORDINATE ][ pntIdx ] = y;
  m_coordinates[ Z_COORDINATE ][ pntIdx ] = z;
}

//------------------------------------------------------------------------------
double MeshCoordinates::getCoordinate( int pntIdx, int dim )
{
  SLIC_ASSERT(  dim < m_ndims );
  SLIC_ASSERT(  ( pntIdx >= 0 ) && ( pntIdx < this->getNumberOfPoints() ) );
  return m_coordinates[ dim ][ pntIdx ];
}

//------------------------------------------------------------------------------
double* MeshCoordinates::getCoordinateArray(int dim)
{
  SLIC_ASSERT( dim < m_ndims );
  return &(m_coordinates[ dim ][ 0 ]);
}

//------------------------------------------------------------------------------
int MeshCoordinates::getNumberOfPoints() const
{
  SLIC_ASSERT( m_ndims >= 1 );
  return m_coordinates[ X_COORDINATE ].size();
}

//------------------------------------------------------------------------------
void MeshCoordinates::initialize()
{
  SLIC_ASSERT( m_ndims >= 1 );

  m_coordinates.resize( m_ndims );
  for ( int i=0 ; i < m_ndims ; ++i )
  {
    m_coordinates[ i ].reserve( 100 );
  }
}

//------------------------------------------------------------------------------
void MeshCoordinates::initialize( int npoints )
{
  SLIC_ASSERT( m_ndims >= 1 );

  m_coordinates.resize( m_ndims );
  for ( int i=0 ; i < m_ndims ; ++i )
  {
    m_coordinates[ i ].resize( npoints );
  }
}

} /* namespace mint */
} /* namespace axom */
