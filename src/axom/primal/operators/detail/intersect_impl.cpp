// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

// Axom includes
#include "axom/primal/operators/detail/intersect_impl.hpp"
#include "axom/slic/interface/slic.hpp"

namespace axom
{
namespace primal
{
namespace detail
{

// -----------------------------------------------------------------------------
bool intersectOnePermutedTriangle(
  const Point3 &p1, const Point3 &q1, const Point3 &r1,
  const Point3 &p2, const Point3 &q2, const Point3 &r2,
  double dp2, double dq2, double dr2,  Vector3 &normal,
  bool includeBoundary, double EPS)
{
  /* Step 4: repeat Step 3, except doing it for triangle 2
     instead of triangle 1 */
  if (isGt(dp2, 0.0, EPS))
  {
    if (isGt(dq2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1,r1,q1,r2,p2,q2,includeBoundary,
                                           EPS);
    }
    else if (isGt(dr2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1,r1,q1,q2,r2,p2,includeBoundary,
                                           EPS);
    }
    else
    {
      return intersectTwoPermutedTriangles(p1,q1,r1,p2,q2,r2,includeBoundary,
                                           EPS);
    }
  }
  else if (isLt(dp2,  0.0, EPS))
  {
    if (isLt(dq2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1,q1,r1,r2,p2,q2,includeBoundary,
                                           EPS);
    }
    else if (isLt(dr2, 0.0, EPS))
    {
      return intersectTwoPermutedTriangles(p1,q1,r1,q2,r2,p2,includeBoundary,
                                           EPS);
    }
    else
    {
      return intersectTwoPermutedTriangles(p1,r1,q1,p2,q2,r2,includeBoundary,
                                           EPS);
    }
  }
  else
  {
    if (isLt(dq2, 0.0, EPS))
    {
      if (isGeq(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1,r1,q1,q2,r2,p2,includeBoundary,
                                             EPS);
      }
      else
      {
        return intersectTwoPermutedTriangles(p1,q1,r1,p2,q2,r2,includeBoundary,
                                             EPS);
      }
    }
    else if (isGt(dq2, 0.0, EPS))
    {
      if (isGt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1,r1,q1,p2,q2,r2,includeBoundary,
                                             EPS);
      }
      else
      {
        return intersectTwoPermutedTriangles(p1,q1,r1,q2,r2,p2,includeBoundary,
                                             EPS);
      }
    }
    else
    {
      if (isGt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1,q1,r1,r2,p2,q2,includeBoundary,
                                             EPS);
      }
      else if (isLt(dr2, 0.0, EPS))
      {
        return intersectTwoPermutedTriangles(p1,r1,q1,r2,p2,q2,includeBoundary,
                                             EPS);
      }
      else
      {
        return intersectCoplanar3DTriangles(p1,q1,r1,p2,q2,r2,normal,
                                            includeBoundary, EPS);
      }
    }
  }
}

// -----------------------------------------------------------------------------
bool intersectCoplanar3DTriangles(const Point3& p1,
                                  const Point3& q1,
                                  const Point3& r1,
                                  const Point3& p2,
                                  const Point3& q2,
                                  const Point3& r2,
                                  Vector3 normal,
                                  bool includeBoundary,
                                  double EPS)
{
  /* Co-planar triangles are projected onto the axis that maximizes their
     area and the 2d intersection used to check if they intersect.
   */

  //find triangle with maximum area:
  for (int i=0 ; i<3 ; i++)
  {
    normal[i] = std::abs(normal[i]);
  }

  if ((isGt(normal[0], normal[2], EPS)) && (isGeq(normal[0], normal[1], EPS)))
  {
    //if x projection area greatest, project on YZ and return 2D checker

    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[2],q1[1]),
                                       Point2::make_point(p1[2],p1[1]),
                                       Point2::make_point(r1[2],r1[1]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[2],q2[1]),
                                       Point2::make_point(p2[2],p2[1]),
                                       Point2::make_point(r2[2],r2[1]));

    return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
  }
  else if (isGt(normal[1],normal[2], EPS) && isGeq(normal[1],normal[0], EPS))
  {
    //if y projection area greatest, project on XZ and return 2D checker
    const Triangle2 t1_2da = Triangle2(Point2::make_point(q1[0],q1[2]),
                                       Point2::make_point(p1[0],p1[2]),
                                       Point2::make_point(r1[0],r1[2]));

    const Triangle2 t2_2da = Triangle2(Point2::make_point(q2[0],q2[2]),
                                       Point2::make_point(p2[0],p2[2]),
                                       Point2::make_point(r2[0],r2[2]));

    return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
  }

  //if z projection area greatest, project on XY and return 2D checker
  const Triangle2 t1_2da = Triangle2(Point2::make_point(p1[0],p1[1]),
                                     Point2::make_point(q1[0],q1[1]),
                                     Point2::make_point(r1[0],r1[1]));

  const Triangle2 t2_2da = Triangle2(Point2::make_point(p2[0],p2[1]),
                                     Point2::make_point(q2[0],q2[1]),
                                     Point2::make_point(r2[0],r2[1]));

  return TriangleIntersection2D(t1_2da, t2_2da, includeBoundary, EPS);
}

// -----------------------------------------------------------------------------
bool TriangleIntersection2D(const Triangle2& t1,
                            const Triangle2& t2,
                            bool includeBoundary,
                            double EPS)
{
  if (isLt(twoDcross(t1[0],t1[1],t1[2]),0.0, EPS))
  {
    if ((isLt(twoDcross(t2[0], t2[1], t2[2]),0.0, EPS)))
    {
      return intersectPermuted2DTriangles(t1[0], t1[2], t1[1],
                                          t2[0], t2[2], t2[1],
                                          includeBoundary, EPS);
    }
    else
    {
      return intersectPermuted2DTriangles(t1[0], t1[2], t1[1],
                                          t2[0], t2[1], t2[2],
                                          includeBoundary, EPS);
    }
  }
  else
  {
    if (isLt(twoDcross(t2[0], t2[1], t2[2]),0.0, EPS))
    {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2],
                                          t2[0], t2[2], t2[1],
                                          includeBoundary, EPS);
    }
    else
    {
      return intersectPermuted2DTriangles(t1[0], t1[1], t1[2],
                                          t2[0], t2[1], t2[2],
                                          includeBoundary, EPS);
    }
  }
}

// -----------------------------------------------------------------------------
bool intersectPermuted2DTriangles(const Point2& p1,
                                  const Point2& q1,
                                  const Point2& r1,
                                  const Point2& p2,
                                  const Point2& q2,
                                  const Point2& r2,
                                  bool includeBoundary,
                                  double EPS)
{
  // Step 2: Orient triangle 2 to be counter clockwise and break the problem
  // into two generic cases (where we test the vertex for intersection or the
  // edges).
  //
  // See paper at https://hal.inria.fr/inria-00072100/document for more details

  if (isGpeq(twoDcross(p2,q2,p1), 0.0, includeBoundary, EPS))
  {
    if (isGpeq(twoDcross(q2,r2,p1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(r2,p2,p1), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        return checkEdge(p1,q1,r1,p2,r2,includeBoundary, EPS); //T1 clockwise
      }
    }
    else
    {
      if (isGpeq(twoDcross(r2,p2,p1), 0.0, includeBoundary, EPS))
      {
        //5 region decomposition with p1 in the +-- region
        return checkEdge(p1,q1,r1,r2,q2,includeBoundary, EPS);
      }
      else
      {
        return checkVertex(p1,q1,r1,p2,q2,r2,includeBoundary, EPS);
      }
    }
  }
  else
  {
    if (isGpeq(twoDcross(q2,r2,p1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(r2,p2,p1), 0.0, includeBoundary, EPS))
      {
        //four region decomposition.  ++- region
        return checkEdge(p1,q1,r1,q2,p2,includeBoundary, EPS);
      }
      else
      {
        return checkVertex(p1,q1,r1,q2,r2,p2,includeBoundary, EPS);
      }
    }
    else
    {
      return checkVertex(p1,q1,r1,r2,p2,q2,includeBoundary, EPS);
    }
  }
}

// -----------------------------------------------------------------------------
bool checkEdge(const Point2& p1,
               const Point2& q1,
               const Point2& r1,
               const Point2& p2,
               const Point2& r2,
               bool includeBoundary,
               double EPS)
{
  if (isGpeq(twoDcross(r2, p2, q1), 0.0, includeBoundary, EPS))
  {
    if (isGpeq(twoDcross(r2, p1, q1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(p1, p2, q1), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        if (isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS) &&
            isGpeq(twoDcross(q1, r1, p2), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      return false;
    }
  }
  else
  {
    if (isGpeq(twoDcross(r2, p2, r1), 0.0, includeBoundary, EPS) &&
        isGpeq(twoDcross(q1, r1, r2), 0.0, includeBoundary, EPS) &&
        isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
}

// -----------------------------------------------------------------------------
inline bool checkVertex(const Point2& p1,
                        const Point2& q1,
                        const Point2& r1,
                        const Point2& p2,
                        const Point2& q2,
                        const Point2& r2,
                        bool includeBoundary,
                        double EPS)
{
  if (isGpeq(twoDcross(r2, p2, q1), 0.0, includeBoundary, EPS))
  {
    if (isGpeq(twoDcross(q2, r2, q1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(p1, p2, q1), 0.0, includeBoundary, EPS))
      {
        if (isLpeq(twoDcross(p1, q2, q1), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
      else
      {
        if (isGpeq(twoDcross(p1, p2, r1), 0.0, includeBoundary, EPS) &&
            isGpeq(twoDcross(r2, p2, r1), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      if (isLpeq(twoDcross(p1, q2, q1), 0.0, includeBoundary, EPS) &&
          isGpeq(twoDcross(q2, r2, r1), 0.0, includeBoundary, EPS) &&
          isGpeq(twoDcross(q1, r1, q2), 0.0, includeBoundary, EPS))
      {
        return true;
      }
      else
      {
        return false;
      }
    }
  }
  else
  {
    if (isGpeq(twoDcross(r2, p2, r1), 0.0, includeBoundary, EPS))
    {
      if (isGpeq(twoDcross(q1, r1, r2), 0.0, includeBoundary, EPS))
      {
        if (isGpeq(twoDcross(r1, p1, p2), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
      else
      {
        if (isGpeq(twoDcross(q1, r1, q2), 0.0, includeBoundary, EPS) &&
            isGpeq(twoDcross(q2, r2, r1), 0.0, includeBoundary, EPS))
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
    else
    {
      return false;
    }
  }
}

// -----------------------------------------------------------------------------
bool intervalsDisjoint(double d0, double d1, double d2, double r)
{
  if (d1 < d0)
  {
    std::swap(d1,d0);  // d0 < d1
  }
  if (d2 > d1)
  {
    std::swap(d2,d1);  // d1 is max(d0,d1,d2)
  }
  else if (d2 < d0)
  {
    std::swap(d2,d0);  // d0 is min(d0,d1,d2)

  }
  SLIC_ASSERT(  d0 <= d1 && d0 <= d2);
  SLIC_ASSERT(  d1 >= d0 && d1 >= d2);

  return d1 < -r || d0 > r;
}

} /* end namespace detail */
} /* end namespace primal */
} /* end namespace axom */
