// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file intersect.hpp
 *
 * \brief Consists of functions to test intersection among geometric primitives.
 */

#ifndef INTERSECTION_BEZIER_HPP_
#define INTERSECTION_BEZIER_HPP_

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/OrientedBoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Ray.hpp"
#include "axom/primal/geometry/Segment.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Triangle.hpp"
#include "axom/primal/geometry/BezierCurve.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/core/utilities/Utilities.hpp"

#include "axom/primal/operators/intersect.hpp"

namespace axom
{
namespace primal
{

/*!
 * \brief Tests if Bezier Curves c1 and c2 intersect.
 * \return status true iff c1 intersects with c2, otherwise false.
 *
 * \param c1, c2 BezierCurve objects to intersect
 * \param sp, tp vector of type T parameter space intersection points (t-values and s-values)
 *  for c1 and c2, respectively
 */
template < typename T, int NDIMS>
bool intersect_bezier( const BezierCurve< T, NDIMS>& c1,
    const BezierCurve< T, NDIMS>& c2,
    std::vector< T >& sp,
    std::vector< T >& tp)
  {
    int ord1=c1.getOrder();
    int ord2=c2.getOrder();
  if ( c1.is_linear(1e-31) && c2.is_linear(1e-31))
  {
    T s;
    T t;
    intersect_2d_linear(c1[0],c1[ord1],c2[0],c2[ord2],s,t);
    if (s>=0.0 && s<1.0 && t >=0.0 && t<1.0)
    {
      sp.push_back(s);
      tp.push_back(t);
    }
  }
  else
  {
    BezierCurve< T, NDIMS> c3(ord1); BezierCurve< T, NDIMS> c4(ord1);
    c1.split_bezier(.5,c3,c4);
    
    std::vector<Point<T, NDIMS>> Ptv3 = c3.getControlPoints();
    BoundingBox < T, NDIMS> b3(Ptv3.data(),ord1+1);
    
    std::vector<Point<T, NDIMS>> Ptv4 = c4.getControlPoints();
    BoundingBox < T, NDIMS> b4(Ptv4.data(),ord1+1);
    
    std::vector<Point<T, NDIMS>> Ptv2 = c2.getControlPoints();
    BoundingBox < T, NDIMS> b2(Ptv2.data(),ord1+1);
    bool checkint = false;
    int intcount = sp.size();
    if (intersect(b3,b2))
    {
      intcount=sp.size();
      checkint=intersect_bezier(c2,c3,tp,sp);
      if (checkint)
      {
        for (int i=intcount; i<static_cast<int>(sp.size()); i++)
        {
        sp[i] = 0.5*sp[i];
        }
      }
    }
    if (intersect(b4,b2))
    {
      intcount=sp.size();
      checkint = intersect_bezier(c2,c4,tp,sp);
      if (checkint)
      {
        for (int i=intcount; i<static_cast<int>(sp.size()); i++)
        {
          sp[i] = .5+0.5*sp[i];
        }
      }
    }
  }

  if (sp.size()>0)
  {
    return true;
  }
  else
  {
    return false;
  }
  }

/*!
 * \brief Intersects two segments defined by their end points (a,d) and (c,b)
 * 
 * \param [in] a,d,c,b the endpoints of the segments
 * \param [out] The parametrized s and t values at which intersection occurs (note these could be outside [0,1]
 */

template < typename T, int NDIMS>
void intersect_2d_linear( const Point<T,NDIMS> &a, const Point<T,NDIMS> &d, const Point<T,NDIMS> &c, const Point<T,NDIMS> &b, T &s, T &t)
{
  T determ = a[1]*b[0]-a[0]*b[1]-a[1]*c[0]+a[0]*c[1]+b[1]*d[0]-c[1]*d[0]-b[0]*d[1]+c[0]*d[1];
  s = (1.0/determ)*((b[1]-c[1])*(b[0]-a[0])+(c[0]-b[0])*(b[1]-a[1]));
  t = 1.0-(1.0/determ)*((a[1]-d[1])*(b[0]-a[0])+(d[0]-a[0])*(b[1]-a[1]));
}
} // namespace primal
} // namespace axom

#endif // PRIMAL_INTERSECTION_BEZIER_HPP_
