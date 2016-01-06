/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file Ray.hpp
 *
 * \date Jan 5, 2016
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#ifndef RAY_HPP_
#define RAY_HPP_

#include "quest/Point.hpp"
#include "quest/Segment.hpp"
#include "quest/Vector.hpp"

#include "slic/slic.hpp"

namespace quest {

/*!
 *******************************************************************************
 * \class
 *
 * \brief Represents a ray, \f$ R(t) \in \mathcal{R}^d \f$ , defined by an
 *  origin point, \f$ P \f$ and a normalized direction vector, \f$ \vec{d} \f$,
 *  \f$ \ni R(t)= P + t\vec{d} \forall t \ge 0 \f$
 *******************************************************************************
 */
template < typename T, int NDIMS >
class Ray
{
public:
    typedef Point< T,NDIMS >   PointType;
    typedef Segment< T,NDIMS > SegmentType;
    typedef Vector< T,NDIMS >  VectorType;

public:

    /*!
     ***************************************************************************
     * \brief Constructs a ray object with the given origin and direction.
     * \param [in] origin the origin of the ray.
     * \param [in] direction the direction of the ray.
     * \pre direction.norm()==1
     ***************************************************************************
     */
    Ray( const PointType& origin, const VectorType& direction );

    /*!
     ***************************************************************************
     * \brief Constructs a ray object from a directed segment.
     * \param [in] S user-supplied segment.
     ***************************************************************************
     */
    explicit Ray( const SegmentType& S );

    /*!
     ***************************************************************************
     * \brief Ray Destructor.
     ***************************************************************************
     */
    ~Ray();

    /*!
     ***************************************************************************
     * \brief Returns the point of origin of this Ray instance.
     * \return origin a point instance corresponding to the origin of the ray.
     ***************************************************************************
     */
    const PointType& origin() const { return m_origin; };

    /*!
     ***************************************************************************
     * \brief Returns a point along the ray at the
     * \param [in] t the value at which to
     * \return
     ***************************************************************************
     */
    PointType at( const T& t) const;

    /*!
     ***************************************************************************
     * \brief Returns the direction vector of this Ray instance.
     * \return direction the direction vector of the ray.
     * \post direction.norm()==1
     ***************************************************************************
     */
    const VectorType& direction() const { return m_direction; };

private:

    /*!
     ***************************************************************************
     * \brief Default Constructor. Does nothing.
     * \note Made private to prevent its use in application code.
     ***************************************************************************
     */
    Ray() { };

    PointType m_origin;
    VectorType m_direction;
};

} /* namespace quest */

//------------------------------------------------------------------------------
// Ray Implementation
//------------------------------------------------------------------------------

namespace quest {

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Ray< T,NDIMS >::Ray( const PointType& origin, const VectorType& direction ) :
    m_origin( origin ),
    m_direction( direction )
{
   SLIC_ASSERT( direction.norm()==1.0 );
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Ray< T,NDIMS >::Ray( const SegmentType& S )
{
   m_origin = S.source();

   VectorType dir( S.source(), S.target() );
   m_direction = dir.unitVector();
}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
Ray< T,NDIMS >::~Ray()
{

}

//------------------------------------------------------------------------------
template < typename T, int NDIMS >
inline Point<T,NDIMS> Ray< T,NDIMS >::at( const T& t ) const
{
   PointType p;
   for ( int i=0; i < NDIMS; ++i ) {
       p[i] = m_origin[i] + t*m_direction[i];
   }
   return ( p );
}

} /* end namespace quest */

#endif /* RAY_HPP_ */
