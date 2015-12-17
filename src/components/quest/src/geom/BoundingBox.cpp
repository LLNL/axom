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


namespace quest
{

//------------------------------------------------------------------------------
template<int DIM>
BoundingBox<DIM>& BoundingBox<DIM>::operator=(const BoundingBox& rhs )
{

  if ( this != &rhs ) {
    m_min = rhs.m_min;
    m_max = rhs.m_max;
  }

  return *this;
}

//------------------------------------------------------------------------------
template<int DIM>
template<typename T>
bool BoundingBox<DIM>::hasPoint(const Point<T,DIM>& otherPt) const
{
    for(int dim = 0; dim < DIM; ++dim)
    {
        if( otherPt[dim] < m_min[dim] || otherPt[dim] >  m_max[dim])
            return false;
    }
    return true;
}


//------------------------------------------------------------------------------
template<int DIM>
bool BoundingBox<DIM>::isValid() const
{
  for(int dim = 0; dim < DIM; ++dim)
  {
      if( m_min[dim] <  m_max[dim])
          return false;
  }
  return true;

}
} /* namespace quest */
