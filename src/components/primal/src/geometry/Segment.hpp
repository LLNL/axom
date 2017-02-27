/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

#ifndef SEGMENT_HPP_
#define SEGMENT_HPP_

#include "primal/Point.hpp"

namespace axom {
namespace primal {

/*!
 *******************************************************************************
 * \class
 *
 * \brief Represents a directed straight line connecting two points
 *  \f$ A,B \in \mathcal{R}^n \f$ where point A is referred to as the source
 *  point and point B is the target.
 *******************************************************************************
 */
template < typename T,int NDIMS >
class Segment
{
public:
  typedef Point< T,NDIMS > PointType;

public:

  /*!
   *****************************************************************************
   * \brief Creates a segment instance from point A to point B.
   * \param A user-supplied source point
   * \param B user-supplied target point
   *****************************************************************************
   */
  Segment(const PointType& A, const PointType& B);

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  ~Segment();

  /*!
   *****************************************************************************
   * \brief Returns the source point of the segment.
   * \return s the source point of the segment.
   *****************************************************************************
   */
  const PointType& source() const { return m_source; };

  /*!
   *****************************************************************************
   * \brief Returns the target point of the segment.
   * \return t the target point of the segment.
   *****************************************************************************
   */
  const PointType& target() const { return m_target; };

private:

  /*!
   *****************************************************************************
   * \brief Default Constructor. Does nothing.
   * \note Made private to prevent its use in application code.
   *****************************************************************************
   */
  Segment() { };

  PointType m_source;
  PointType m_target;
};

} /* namespace primal */

} /* namespace axom */

//------------------------------------------------------------------------------
//  Segment Implementation
//------------------------------------------------------------------------------
namespace axom {
namespace primal {

template < typename T, int NDIMS >
Segment< T,NDIMS >::Segment(const PointType& A, const PointType& B):
  m_source( A ),
  m_target( B )
{}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Segment< T,NDIMS >::~Segment() { }

} /* namespacce primal */

} /* namespace axom */

#endif /* SEGMENT_HPP_ */
