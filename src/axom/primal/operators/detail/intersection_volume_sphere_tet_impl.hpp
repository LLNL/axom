// Copyright (c) Lawrence Livermore National Security, LLC and other
// Axom Project Contributors. See top-level LICENSE and COPYRIGHT
// files for dates and other details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#pragma once

/*!
 * \file intersection_volume_sphere_tet_impl.hpp
 *
 * \brief Implements exact sphere-tetrahedron overlap volume.
 *
 * The algorithm follows S. Strobl, A. Formella, and T. Poeschel,
 * "Exact calculation of the overlap volume of spheres and mesh elements,"
 * Journal of Computational Physics 311 (2016), 158-172.
 * https://doi.org/10.1016/j.jcp.2016.02.003
 *
 * The implementation adapts robustness details from the authors' MIT-licensed
 * reference at https://github.com/severinstrobl/overlap, revision
 * 7c5db4ed03653b560eafe0914397205c96e35cc2.
 *
 * MIT License
 *
 * Copyright (c) 2015-2025 Severin Strobl
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "axom/core/NumericLimits.hpp"
#include "axom/core/StackArray.hpp"
#include "axom/primal/geometry/Sphere.hpp"
#include "axom/primal/geometry/Tetrahedron.hpp"

#include <cmath>

namespace axom
{
namespace primal
{
namespace detail
{
namespace sphere_tet
{

using Vector3 = primal::Vector<double, 3>;

constexpr int NUM_VERTICES = 4;
constexpr int NUM_EDGES = 6;
constexpr int NUM_FACES = 4;
constexpr double TINY_EPSILON = 4.0 * axom::numeric_limits<double>::epsilon();
constexpr double LARGE_EPSILON = 1.e-10;

/*!
 * \brief Returns an endpoint vertex for a tetrahedron edge.
 * \param [in] edge the edge index
 * \param [in] endpoint the local endpoint index
 * \return the tetrahedron vertex index
 */
AXOM_HOST_DEVICE inline int edgeVertex(int edge, int endpoint)
{
  if(edge < 3)
  {
    return endpoint == 0 ? edge : (edge + 1) % 3;
  }
  return endpoint == 0 ? edge - 3 : 3;
}

/*!
 * \brief Returns an incident face for a tetrahedron edge.
 * \param [in] edge the edge index
 * \param [in] side the local face index
 * \return the tetrahedron face index
 */
AXOM_HOST_DEVICE inline int edgeFace(int edge, int side)
{
  if(side == 0)
  {
    return edge < 3 ? 0 : (edge < 5 ? 1 : 2);
  }

  switch(edge)
  {
  case 0:
    return 1;
  case 1:
  case 4:
    return 2;
  default:
    return 3;
  }
}

/*!
 * \brief Returns an incident edge for a tetrahedron vertex.
 * \param [in] vertex the vertex index
 * \param [in] localEdge the local incident-edge index
 * \return the tetrahedron edge index
 */
AXOM_HOST_DEVICE inline int vertexEdge(int vertex, int localEdge)
{
  switch(vertex)
  {
  case 0:
    return localEdge == 0 ? 0 : (localEdge == 1 ? 2 : 3);
  case 1:
    return localEdge == 0 ? 0 : (localEdge == 1 ? 1 : 4);
  case 2:
    return localEdge == 0 ? 1 : (localEdge == 1 ? 2 : 5);
  default:
    return localEdge + 3;
  }
}

/*!
 * \brief Returns the endpoint index of a vertex on an incident edge.
 * \param [in] vertex the vertex index
 * \param [in] localEdge the local incident-edge index
 * \return the local endpoint index
 */
AXOM_HOST_DEVICE inline int vertexEdgeEnd(int vertex, int localEdge)
{
  if(vertex == 0)
  {
    return localEdge == 1 ? 1 : 0;
  }
  return vertex == 3 || localEdge == 0 ? 1 : 0;
}

/*!
 * \brief Returns an incident face for a tetrahedron vertex.
 * \param [in] vertex the vertex index
 * \param [in] localFace the local incident-face index
 * \return the tetrahedron face index
 */
AXOM_HOST_DEVICE inline int vertexFace(int vertex, int localFace)
{
  switch(vertex)
  {
  case 0:
    return localFace < 2 ? localFace : 3;
  case 1:
    return localFace;
  case 2:
    return localFace == 0 ? 0 : localFace + 1;
  default:
    return localFace == 0 ? 1 : (localFace == 1 ? 3 : 2);
  }
}

/*!
 * \brief Returns an edge index in the local vertex-correction triangle.
 * \param [in] localFace the local face index
 * \param [in] side the local edge endpoint index
 * \return the vertex-correction point index
 */
AXOM_HOST_DEVICE inline int vertexFaceEdge(int localFace, int side)
{
  return side == 0 ? (localFace == 2 ? 1 : 0) : (localFace == 0 ? 1 : 2);
}

/*!
 * \brief Returns a vertex for an oriented tetrahedron face.
 * \param [in] face the face index
 * \param [in] localVertex the local face-vertex index
 * \return the tetrahedron vertex index
 */
AXOM_HOST_DEVICE inline int faceVertex(int face, int localVertex)
{
  switch(face)
  {
  case 0:
    return localVertex == 0 ? 2 : (localVertex == 1 ? 1 : 0);
  case 1:
    return localVertex == 2 ? 3 : localVertex;
  case 2:
    return localVertex == 0 ? 1 : (localVertex == 1 ? 2 : 3);
  default:
    return localVertex == 0 ? 2 : (localVertex == 1 ? 0 : 3);
  }
}

/*!
 * \brief Computes the dot product of two vectors.
 * \param [in] lhs the first vector
 * \param [in] rhs the second vector
 * \return the dot product
 */
AXOM_HOST_DEVICE inline double dot(const Vector3& lhs, const Vector3& rhs)
{
  return Vector3::dot_product(lhs, rhs);
}

/*!
 * \brief Computes the cross product of two vectors.
 * \param [in] lhs the first vector
 * \param [in] rhs the second vector
 * \return the cross product
 */
AXOM_HOST_DEVICE inline Vector3 cross(const Vector3& lhs, const Vector3& rhs)
{
  return Vector3::cross_product(lhs, rhs);
}

/*!
 * \brief Computes a vector norm without avoidable overflow or underflow.
 * \param [in] vector the vector
 * \return the Euclidean norm
 */
AXOM_HOST_DEVICE inline double stableNorm(const Vector3& vector)
{
  double scale = 0.0;
  for(int dim = 0; dim < 3; ++dim)
  {
    scale = axom::utilities::max(scale, ::fabs(vector[dim]));
  }
  if(scale == 0.0)
  {
    return 0.0;
  }

  double sum = 0.0;
  for(int dim = 0; dim < 3; ++dim)
  {
    const double value = vector[dim] / scale;
    sum += value * value;
  }
  return scale * ::sqrt(sum);
}

/*!
 * \brief Normalizes a vector using its stable norm.
 * \param [in] vector the vector
 * \return the normalized vector, or zero when \p vector is zero
 */
AXOM_HOST_DEVICE inline Vector3 stableNormalized(const Vector3& vector)
{
  const double norm = stableNorm(vector);
  return norm > 0.0 ? vector / norm : Vector3 {};
}

/*! \brief Stores a double-double value for robust orientation predicates. */
class DoubleDouble
{
public:
  AXOM_HOST_DEVICE explicit DoubleDouble(double value = 0.0) : m_high(value), m_low(0.0) { }

  AXOM_HOST_DEVICE DoubleDouble(double high, double low) : m_high(high), m_low(low) { }

  AXOM_HOST_DEVICE friend DoubleDouble operator+(const DoubleDouble& lhs, const DoubleDouble& rhs)
  {
    const double sum = lhs.m_high + rhs.m_high;
    const double value = sum - lhs.m_high;
    const double error = (lhs.m_high - (sum - value)) + (rhs.m_high - value) + lhs.m_low + rhs.m_low;
    const double result = sum + error;
    return DoubleDouble(result, (sum - result) + error);
  }

  AXOM_HOST_DEVICE friend DoubleDouble operator-(const DoubleDouble& lhs, const DoubleDouble& rhs)
  {
    return lhs + DoubleDouble(-rhs.m_high, -rhs.m_low);
  }

  AXOM_HOST_DEVICE friend DoubleDouble operator*(const DoubleDouble& lhs, const DoubleDouble& rhs)
  {
    const double product = lhs.m_high * rhs.m_high;
    const double error = axom::utilities::fma(lhs.m_high, rhs.m_high, -product);
    const double correction =
      error + lhs.m_high * rhs.m_low + lhs.m_low * rhs.m_high + lhs.m_low * rhs.m_low;
    const double result = product + correction;
    return DoubleDouble(result, (product - result) + correction);
  }

  AXOM_HOST_DEVICE double value() const { return m_high + m_low; }

private:
  double m_high;
  double m_low;
};

/*!
 * \brief Evaluates the two-dimensional orientation determinant with extended precision.
 * \param [in] ax x coordinate of the first point
 * \param [in] ay y coordinate of the first point
 * \param [in] bx x coordinate of the second point
 * \param [in] by y coordinate of the second point
 * \param [in] cx x coordinate of the third point
 * \param [in] cy y coordinate of the third point
 * \return the signed orientation of the three input points
 */
AXOM_HOST_DEVICE inline double orient2D(double ax, double ay, double bx, double by, double cx, double cy)
{
  const DoubleDouble x0(ax);
  const DoubleDouble y0(ay);
  const DoubleDouble x1(bx);
  const DoubleDouble y1(by);
  const DoubleDouble x2(cx);
  const DoubleDouble y2(cy);
  return ((x0 - x2) * (y1 - y2) - (y0 - y2) * (x1 - x2)).value();
}

/*!
 * \brief Computes a robust unit normal for a triangle.
 * \param [in] a the first triangle point
 * \param [in] b the second triangle point
 * \param [in] c the third triangle point
 * \return the unit normal
 */
AXOM_HOST_DEVICE inline Vector3 robustTriangleNormal(const Vector3& a,
                                                     const Vector3& b,
                                                     const Vector3& c)
{
  const double xy = orient2D(a[0], a[1], b[0], b[1], c[0], c[1]);
  const double yz = orient2D(a[1], a[2], b[1], b[2], c[1], c[2]);
  const double zx = orient2D(a[2], a[0], b[2], b[0], c[2], c[0]);
  return stableNormalized(Vector3 {yz, zx, xy});
}

/*!
 * \brief Computes the angle between normalized vectors robustly near zero and pi.
 * \param [in] lhs the first normalized vector
 * \param [in] rhs the second normalized vector
 * \return the angle in radians
 */
AXOM_HOST_DEVICE inline double angle(const Vector3& lhs, const Vector3& rhs)
{
  if(dot(lhs, rhs) < 0.0)
  {
    return M_PI - 2.0 * ::asin(0.5 * stableNorm(-rhs - lhs));
  }
  return 2.0 * ::asin(0.5 * stableNorm(rhs - lhs));
}

struct Face
{
  Vector3 center;
  Vector3 normal;
};

struct Tet
{
  axom::StackArray<Vector3, NUM_VERTICES> vertices;
  axom::StackArray<Face, NUM_FACES> faces;
  Vector3 center;
  double volume;
};

/*!
 * \brief Constructs an oriented triangular face.
 * \param [in] a the first face vertex
 * \param [in] b the second face vertex
 * \param [in] c the third face vertex
 * \return the face with its center and unit normal initialized
 */
AXOM_HOST_DEVICE inline Face makeFace(const Vector3& a, const Vector3& b, const Vector3& c)
{
  Face face;
  face.center = (a + b + c) / 3.0;
  face.normal = stableNormalized(cross(b - a, c - a));
  return face;
}

/*!
 * \brief Constructs a positively oriented tetrahedron and its outward faces.
 * \param [in] vertices the tetrahedron vertices
 * \return the initialized tetrahedron
 */
AXOM_HOST_DEVICE inline Tet makeTet(axom::StackArray<Vector3, NUM_VERTICES> vertices)
{
  if(Vector3::scalar_triple_product(vertices[1] - vertices[0],
                                    vertices[2] - vertices[0],
                                    vertices[3] - vertices[0]) < 0.0)
  {
    axom::utilities::swap(vertices[1], vertices[2]);
  }

  Tet tet;
  tet.vertices = vertices;
  tet.faces[0] = makeFace(vertices[2], vertices[1], vertices[0]);
  tet.faces[1] = makeFace(vertices[0], vertices[1], vertices[3]);
  tet.faces[2] = makeFace(vertices[1], vertices[2], vertices[3]);
  tet.faces[3] = makeFace(vertices[2], vertices[0], vertices[3]);
  tet.center = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / 4.0;
  tet.volume = ::fabs(Vector3::scalar_triple_product(vertices[0] - vertices[3],
                                                     vertices[1] - vertices[3],
                                                     vertices[2] - vertices[3])) /
    6.0;
  return tet;
}

/*!
 * \brief Determines whether a point projects inside a triangular face.
 * \param [in] tet the tetrahedron
 * \param [in] faceIndex the face index
 * \param [in] point the point
 * \return true when the projected point is inside the face
 */
AXOM_HOST_DEVICE inline bool contains(const Tet& tet, int faceIndex, const Vector3& point)
{
  const Face& face = tet.faces[faceIndex];
  const Vector3 projected = point - dot(face.normal, point - face.center) * face.normal;
  for(int i = 0; i < 3; ++i)
  {
    const Vector3& a = tet.vertices[faceVertex(faceIndex, i)];
    const Vector3& b = tet.vertices[faceVertex(faceIndex, (i + 1) % 3)];
    const Vector3 direction = cross(b - a, face.normal);
    if(dot(direction, projected - 0.5 * (a + b)) > 0.0)
    {
      return false;
    }
  }
  return true;
}

/*!
 * \brief Determines whether a tetrahedron contains a point.
 * \param [in] tet the tetrahedron
 * \param [in] point the point
 * \return true when the point is inside every face
 */
AXOM_HOST_DEVICE inline bool contains(const Tet& tet, const Vector3& point)
{
  for(int face = 0; face < NUM_FACES; ++face)
  {
    if(dot(tet.faces[face].normal, point - tet.faces[face].center) > 0.0)
    {
      return false;
    }
  }
  return true;
}

/*!
 * \brief Computes the volume of a unit-sphere cap.
 * \param [in] height the cap height
 * \return the cap volume
 */
AXOM_HOST_DEVICE inline double capVolume(double height)
{
  if(height <= 0.0)
  {
    return 0.0;
  }
  if(height >= 2.0)
  {
    return 4.0 * M_PI / 3.0;
  }
  return (M_PI / 3.0) * height * height * (3.0 - height);
}

struct EdgeRoots
{
  bool intersects {false};
  double first {0.0};
  double second {0.0};
};

/*!
 * \brief Intersects a line with the unit sphere.
 * \param [in] base the line origin
 * \param [in] direction the line direction
 * \return the ordered line parameters when two proper intersections exist
 */
AXOM_HOST_DEVICE inline EdgeRoots intersectUnitSphere(const Vector3& base, const Vector3& direction)
{
  const double a = dot(direction, direction);
  if(a == 0.0)
  {
    return {};
  }

  const double b = 2.0 * dot(direction, base);
  const double c = dot(base, base) - 1.0;
  const double discriminant = b * b - 4.0 * a * c;
  if(discriminant <= 0.0)
  {
    return {};
  }

  const double q = -0.5 * (b + ::copysign(::sqrt(discriminant), b));
  double first = q / a;
  double second = c / q;
  if(first > second)
  {
    axom::utilities::swap(first, second);
  }

  const Vector3 midpoint = base + (0.5 * (first + second)) * direction;
  if(dot(midpoint, midpoint) >= 1.0)
  {
    return {};
  }
  return {true, first, second};
}

/*!
 * \brief Computes the principal unit-sphere wedge term.
 * \param [in] distance distance from the sphere center to the wedge edge
 * \param [in] alpha the wedge angle
 * \return the principal wedge volume
 */
AXOM_HOST_DEVICE inline double principalWedge(double distance, double alpha)
{
  const double sinAlpha = ::sin(alpha);
  const double cosAlpha = ::cos(alpha);
  const double a = distance * sinAlpha;
  const double b = ::sqrt(::fabs(1.0 - distance * distance));
  const double c = distance * cosAlpha;
  return (a * b * c) / 3.0 + a * (a * a / 3.0 - 1.0) * ::atan2(b, c) +
    (2.0 / 3.0) * ::atan2(sinAlpha * b, cosAlpha);
}

/*!
 * \brief Computes a regularized unit-sphere wedge volume.
 * \param [in] distance distance from the sphere center to the wedge edge
 * \param [in] alpha the wedge angle
 * \param [in] axialDistance signed distance along the wedge axis
 * \return the wedge volume
 */
AXOM_HOST_DEVICE inline double regularizedWedge(double distance, double alpha, double axialDistance)
{
  const bool obtuse = alpha > 0.5 * M_PI;
  const double acuteAlpha = obtuse ? M_PI - alpha : alpha;
  const double principal = principalWedge(distance, acuteAlpha);

  if(axialDistance >= 0.0)
  {
    if(obtuse)
    {
      const double height = 1.0 - axialDistance;
      return (M_PI / 3.0) * height * height * (3.0 - height) - principal;
    }
    return principal;
  }

  const double hemisphereVolume = 2.0 * M_PI / 3.0;
  if(obtuse)
  {
    return hemisphereVolume - principal;
  }
  return hemisphereVolume - (capVolume(1.0 + axialDistance) - principal);
}

struct Plane
{
  Vector3 center;
  Vector3 normal;
};

/*!
 * \brief Computes the unit-sphere wedge bounded by two planes.
 * \param [in] first the first bounding plane
 * \param [in] second the second bounding plane
 * \param [in] offset the vector from the sphere center to the planes' intersection line
 * \return the wedge volume
 */
AXOM_HOST_DEVICE inline double generalWedge(const Plane& first,
                                            const Plane& second,
                                            const Vector3& offset)
{
  const double distance = stableNorm(offset);
  if(distance < TINY_EPSILON)
  {
    return (2.0 / 3.0) * (M_PI - angle(first.normal, second.normal));
  }
  if(distance >= 1.0)
  {
    return 0.0;
  }

  const double firstDistance = dot(offset, first.normal);
  const double secondDistance = dot(offset, second.normal);
  if(::fabs(firstDistance) < TINY_EPSILON || ::fabs(secondDistance) < TINY_EPSILON)
  {
    const double alpha = M_PI - angle(first.normal, second.normal);
    const double axial =
      ::fabs(firstDistance) > ::fabs(secondDistance) ? firstDistance : secondDistance;
    return regularizedWedge(distance, alpha, axial);
  }

  Vector3 direction = offset / distance;
  if(distance < LARGE_EPSILON)
  {
    const Vector3 normalCross = stableNormalized(cross(first.normal, second.normal));
    direction = stableNormalized(direction - dot(direction, normalCross) * normalCross);
  }

  double firstAngle = angle(first.normal, direction);
  double secondAngle = angle(second.normal, direction);
  const double firstDirection = dot(direction, offset - first.center);
  const double secondDirection = dot(direction, offset - second.center);

  if(firstDistance >= 0.0 && secondDistance >= 0.0)
  {
    firstAngle = 0.5 * M_PI - ::copysign(firstAngle, firstDirection);
    secondAngle = 0.5 * M_PI - ::copysign(secondAngle, secondDirection);
    return regularizedWedge(distance, firstAngle, firstDistance) +
      regularizedWedge(distance, secondAngle, secondDistance);
  }

  if(firstDistance < 0.0 && secondDistance < 0.0)
  {
    firstAngle = 0.5 * M_PI + ::copysign(1.0, firstDirection) * (firstAngle - M_PI);
    secondAngle = 0.5 * M_PI + ::copysign(1.0, secondDirection) * (secondAngle - M_PI);
    return 4.0 * M_PI / 3.0 -
      (regularizedWedge(distance, firstAngle, -firstDistance) +
       regularizedWedge(distance, secondAngle, -secondDistance));
  }

  firstAngle = 0.5 * M_PI -
    ::copysign(1.0, firstDirection * firstDistance) *
      (firstAngle - (firstDistance < 0.0 ? M_PI : 0.0));
  secondAngle = 0.5 * M_PI -
    ::copysign(1.0, secondDirection * secondDistance) *
      (secondAngle - (secondDistance < 0.0 ? M_PI : 0.0));

  const double firstVolume = regularizedWedge(distance, firstAngle, ::fabs(firstDistance));
  const double secondVolume = regularizedWedge(distance, secondAngle, ::fabs(secondDistance));
  return ::fabs(firstVolume - secondVolume);
}

struct IntersectionData
{
  unsigned int vertices;
  unsigned int edges;
  unsigned int faces;
  axom::StackArray<axom::StackArray<Vector3, 2>, NUM_EDGES> edgePoints;
};

/*!
 * \brief Determines whether an entity is present in a compact bit mask.
 * \param [in] entities the entity bit mask
 * \param [in] index the entity index
 * \return true when the entity is marked
 */
AXOM_HOST_DEVICE inline bool isMarked(unsigned int entities, int index)
{
  return (entities & (1U << index)) != 0;
}

/*!
 * \brief Determines whether a triangular face intersects the unit sphere.
 * \param [in] tet the tetrahedron
 * \param [in] face the face index
 * \return true when the sphere crosses the face interior
 */
AXOM_HOST_DEVICE inline bool faceIntersectsUnitSphere(const Tet& tet, int face)
{
  const double distance = dot(tet.faces[face].normal, -tet.faces[face].center);
  return distance * distance < 1.0 && contains(tet, face, Vector3 {});
}

/*!
 * \brief Finds sphere intersections with tetrahedron vertices, edges, and faces.
 * \param [in] tet the tetrahedron in unit-sphere coordinates
 * \return the marked entities and edge intersection points
 */
AXOM_HOST_DEVICE inline IntersectionData findIntersections(const Tet& tet)
{
  IntersectionData data;
  data.vertices = 0;
  data.edges = 0;
  data.faces = 0;

  for(int edge = 0; edge < NUM_EDGES; ++edge)
  {
    const int firstVertex = edgeVertex(edge, 0);
    const int secondVertex = edgeVertex(edge, 1);
    const Vector3& base = tet.vertices[firstVertex];
    const Vector3 direction = tet.vertices[secondVertex] - base;
    const EdgeRoots roots = intersectUnitSphere(base, direction);
    if(!roots.intersects || roots.first >= 1.0 || roots.second <= 0.0)
    {
      continue;
    }

    data.vertices |= (roots.first < 0.0 ? 1U : 0U) << firstVertex;
    data.vertices |= (roots.second > 1.0 ? 1U : 0U) << secondVertex;
    data.edgePoints[edge][0] = roots.first * direction;
    data.edgePoints[edge][1] = (roots.second - 1.0) * direction;
    data.edges |= 1U << edge;
    data.faces |= 1U << edgeFace(edge, 0);
    data.faces |= 1U << edgeFace(edge, 1);
  }

  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    if(isMarked(data.vertices, vertex))
    {
      const bool allEdgesIntersect = isMarked(data.edges, vertexEdge(vertex, 0)) &&
        isMarked(data.edges, vertexEdge(vertex, 1)) && isMarked(data.edges, vertexEdge(vertex, 2));
      if(!allEdgesIntersect)
      {
        data.vertices &= ~(1U << vertex);
      }
    }
  }

  for(int face = 0; face < NUM_FACES; ++face)
  {
    data.faces |= (faceIntersectsUnitSphere(tet, face) ? 1U : 0U) << face;
  }
  return data;
}

/*!
 * \brief Computes the spherical wedge contribution for an intersected edge.
 * \param [in] tet the tetrahedron
 * \param [in] data the cached sphere intersections
 * \param [in] edge the edge index
 * \return the edge wedge volume
 */
AXOM_HOST_DEVICE inline double edgeWedge(const Tet& tet, const IntersectionData& data, int edge)
{
  const Face& firstFace = tet.faces[edgeFace(edge, 0)];
  const Face& secondFace = tet.faces[edgeFace(edge, 1)];
  const Vector3 midpoint = 0.5 *
    ((data.edgePoints[edge][0] + tet.vertices[edgeVertex(edge, 0)]) +
     (data.edgePoints[edge][1] + tet.vertices[edgeVertex(edge, 1)]));
  return generalWedge({firstFace.center, firstFace.normal},
                      {secondFace.center, secondFace.normal},
                      midpoint);
}

/*!
 * \brief Computes the spherical-cone correction at an intersected vertex.
 * \param [in] tet the tetrahedron
 * \param [in] data the cached sphere intersections
 * \param [in] vertex the vertex index
 * \return the vertex correction volume
 */
AXOM_HOST_DEVICE inline double vertexCorrection(const Tet& tet, const IntersectionData& data, int vertex)
{
  struct IndexedDistance
  {
    int index;
    double squaredDistance;
  };

  axom::StackArray<Vector3, 3> relativePoints;
  axom::StackArray<Vector3, 3> points;
  axom::StackArray<IndexedDistance, 3> distances;
  for(int localEdge = 0; localEdge < 3; ++localEdge)
  {
    const int edge = vertexEdge(vertex, localEdge);
    relativePoints[localEdge] = data.edgePoints[edge][vertexEdgeEnd(vertex, localEdge)];
    points[localEdge] = relativePoints[localEdge] + tet.vertices[vertex];
    distances[localEdge] = {localEdge, dot(relativePoints[localEdge], relativePoints[localEdge])};
  }
  if(distances[1].squaredDistance < distances[0].squaredDistance)
  {
    axom::utilities::swap(distances[0], distances[1]);
  }
  if(distances[2].squaredDistance < distances[1].squaredDistance)
  {
    axom::utilities::swap(distances[1], distances[2]);
  }
  if(distances[1].squaredDistance < distances[0].squaredDistance)
  {
    axom::utilities::swap(distances[0], distances[1]);
  }

  if(distances[1].squaredDistance < distances[2].squaredDistance * LARGE_EPSILON)
  {
    return edgeWedge(tet, data, vertexEdge(vertex, distances[2].index));
  }

  const Vector3 center = (points[0] + points[1] + points[2]) / 3.0;
  Vector3 normal = robustTriangleNormal(relativePoints[0], relativePoints[1], relativePoints[2]);
  if(dot(normal, tet.center - center) > 0.0)
  {
    normal = -normal;
  }

  const Plane conePlane {center, normal};
  double segmentCorrection = 0.0;
  for(int localFace = 0; localFace < 3; ++localFace)
  {
    const Face& face = tet.faces[vertexFace(vertex, localFace)];
    const Vector3 midpoint =
      0.5 * (points[vertexFaceEdge(localFace, 0)] + points[vertexFaceEdge(localFace, 1)]);
    segmentCorrection += generalWedge(conePlane, {face.center, -face.normal}, midpoint);
  }

  const double tipVolume =
    ::fabs(dot(-relativePoints[2],
               cross(relativePoints[0] - relativePoints[2], relativePoints[1] - relativePoints[2]))) /
    6.0;
  const double cap = capVolume(1.0 + dot(normal, -center));
  if(cap < TINY_EPSILON)
  {
    return tipVolume;
  }
  return axom::utilities::max(0.0, tipVolume + cap - segmentCorrection);
}

/*!
 * \brief Computes unit-sphere overlap with a normalized tetrahedron.
 * \param [in] tet the tetrahedron in unit-sphere coordinates
 * \return the overlap volume
 */
AXOM_HOST_DEVICE inline double unitSphereTetOverlap(const Tet& tet)
{
  if(tet.volume == 0.0)
  {
    return 0.0;
  }

  bool allVerticesInside = true;
  bool boxesOverlap = true;
  for(int dim = 0; dim < 3; ++dim)
  {
    double lower = tet.vertices[0][dim];
    double upper = lower;
    for(int vertex = 1; vertex < NUM_VERTICES; ++vertex)
    {
      lower = axom::utilities::min(lower, tet.vertices[vertex][dim]);
      upper = axom::utilities::max(upper, tet.vertices[vertex][dim]);
    }
    boxesOverlap = boxesOverlap && lower <= 1.0 && upper >= -1.0;
  }
  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    allVerticesInside = allVerticesInside && dot(tet.vertices[vertex], tet.vertices[vertex]) <= 1.0;
  }
  if(!boxesOverlap)
  {
    return 0.0;
  }
  if(allVerticesInside)
  {
    return tet.volume;
  }

  const IntersectionData intersections = findIntersections(tet);
  if(intersections.faces == 0 && contains(tet, Vector3 {}))
  {
    return 4.0 * M_PI / 3.0;
  }
  if(intersections.vertices == 0 && intersections.edges == 0 && intersections.faces == 0)
  {
    return 0.0;
  }

  double result = 4.0 * M_PI / 3.0;
  for(int face = 0; face < NUM_FACES; ++face)
  {
    if(isMarked(intersections.faces, face))
    {
      result -= capVolume(1.0 + dot(tet.faces[face].normal, -tet.faces[face].center));
    }
  }
  for(int edge = 0; edge < NUM_EDGES; ++edge)
  {
    if(isMarked(intersections.edges, edge))
    {
      result += edgeWedge(tet, intersections, edge);
    }
  }
  for(int vertex = 0; vertex < NUM_VERTICES; ++vertex)
  {
    if(isMarked(intersections.vertices, vertex))
    {
      result -= vertexCorrection(tet, intersections, vertex);
    }
  }

  const double maximum = axom::utilities::min(4.0 * M_PI / 3.0, tet.volume);
  const double tolerance = ::sqrt(axom::numeric_limits<double>::epsilon()) * maximum;
  if(result < 0.0 && result > -tolerance)
  {
    result = 0.0;
  }
  if(result > maximum && result - maximum < tolerance)
  {
    result = maximum;
  }
  SLIC_ASSERT(result >= 0.0 && result <= maximum);
  return result;
}

}  // namespace sphere_tet
}  // namespace detail

/*!
 * \brief Finds the exact sphere-tetrahedron intersection volume.
 *
 * The computation normalizes the tetrahedron against a unit sphere and uses
 * inclusion-exclusion over spherical caps, edge wedges, and vertex cones.
 *
 * \tparam T the coordinate type
 * \param [in] sphere the sphere
 * \param [in] tet the tetrahedron
 * \return the unsigned intersection volume
 * \note This overload is available in host and device code.
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Sphere<T, 3>& sphere, const Tetrahedron<T, 3>& tet)
{
  const double radius = static_cast<double>(sphere.getRadius());
  if(radius <= 0.0)
  {
    return T {};
  }

  const auto& center = sphere.getCenter();
  axom::StackArray<detail::sphere_tet::Vector3, detail::sphere_tet::NUM_VERTICES> vertices;
  for(int vertex = 0; vertex < detail::sphere_tet::NUM_VERTICES; ++vertex)
  {
    vertices[vertex] = detail::sphere_tet::Vector3 {
      (static_cast<double>(tet[vertex][0]) - static_cast<double>(center[0])) / radius,
      (static_cast<double>(tet[vertex][1]) - static_cast<double>(center[1])) / radius,
      (static_cast<double>(tet[vertex][2]) - static_cast<double>(center[2])) / radius};
  }
  const auto normalizedTet = detail::sphere_tet::makeTet(vertices);
  const double normalizedVolume = detail::sphere_tet::unitSphereTetOverlap(normalizedTet);
  return static_cast<T>(normalizedVolume * radius * radius * radius);
}

/*!
 * \brief Finds the exact tetrahedron-sphere intersection volume.
 *
 * \tparam T the coordinate type
 * \param [in] tet the tetrahedron
 * \param [in] sphere the sphere
 * \return the unsigned intersection volume
 * \note This overload is available in host and device code.
 */
template <typename T>
AXOM_HOST_DEVICE T intersection_volume(const Tetrahedron<T, 3>& tet, const Sphere<T, 3>& sphere)
{
  return intersection_volume(sphere, tet);
}

}  // namespace primal
}  // namespace axom
