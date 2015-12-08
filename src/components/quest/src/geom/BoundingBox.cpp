/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file BoundingBox.cpp
 *
 * \date Dec 9, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "BoundingBox.hpp"

// C/C++ includes
#include <algorithm> // for fill

namespace quest
{

BoundingBox::BoundingBox()
{
  std::fill( m_min, m_min+3, 0.0 );
  std::fill( m_max, m_max+3, 1.0 );
}

//------------------------------------------------------------------------------
BoundingBox::~BoundingBox()
{

}

//------------------------------------------------------------------------------
BoundingBox& BoundingBox::operator=(const BoundingBox& rhs )
{

  if ( this != &rhs ) {
    memcpy( m_min, rhs.m_min, 3*sizeof(double) );
    memcpy( m_max, rhs.m_max, 3*sizeof(double) );
  }

  return *this;
}

//------------------------------------------------------------------------------
bool BoundingBox::hasPoint( double x, double y, double z ) const
{
  bool status = false;

  if ( x >= m_min[0] && x <= m_max[0] &&
       y >= m_min[1] && y <= m_max[1] &&
       z >= m_min[2] && z <= m_max[2]       ) {

       status = true;

  }

  return status;
}

} /* namespace quest */
