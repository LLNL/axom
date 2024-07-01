// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/primal.hpp"
#include "axom/slic.hpp"

#include "axom/core/execution/execution_space.hpp"

#include <math.h>
#include "gtest/gtest.h"

namespace
{
constexpr double EPS = 1e-8;
}

//------------------------------------------------------------------------------
TEST(primal_polygon, empty)
{
  using PolygonType = axom::primal::Polygon<double, 3>;
  PolygonType poly;
  EXPECT_FALSE(poly.isValid());
}

//------------------------------------------------------------------------------
TEST(primal_polygon, winding_number)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Test the specific winding numbers
  EXPECT_EQ(winding_number(PointType {0.25, 0.5}, poly), 1);
  EXPECT_EQ(winding_number(PointType {0.75, 0.5}, poly), -1);
  EXPECT_EQ(winding_number(PointType {0.5, 0.25}, poly), 0);
  EXPECT_EQ(winding_number(PointType {0.5, 0.75}, poly), 0);

  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);
  EXPECT_EQ(winding_number(PointType {-0.1, 0.0}, poly), -2);
  EXPECT_EQ(winding_number(PointType {0.1, 0.0}, poly), -1);
  EXPECT_EQ(winding_number(PointType {-2.0, 0.0}, poly), 0);
  EXPECT_EQ(winding_number(PointType {2.5, 0.0}, poly), 0);

  // Current policy is to return 1 on a boundary when it is included,
  //  0 on boundaries when it is not, as 0 always indicates "exterior"
  const bool includeBoundary = true;
  EXPECT_EQ(winding_number(PointType {0.0, 0.0}, poly, includeBoundary), 1);
  EXPECT_EQ(winding_number(PointType {0.0, 0.0}, poly, !includeBoundary), 0);
}

//------------------------------------------------------------------------------
TEST(primal_polygon, containment)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly(vertices);

  // Simple cases on non-convex polygon
  EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
  EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));

  // Edge cases, where vertex is aligned with edge
  EXPECT_FALSE(in_polygon(PointType {-0.25, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.25, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {-0.25, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.25, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.0, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.0, -0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {-1.0, 0.0}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.0, 1.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {1.0, -0.25}, poly));

  // Test different in/out protocols for polygon with extra loop
  vertices = axom::Array<PointType>({PointType {0, 1},
                                     PointType {0, -1},
                                     PointType {-2, 2},
                                     PointType {2, 2},
                                     PointType {2, -2},
                                     PointType {-2, -2}});

  poly = PolygonType(vertices);

  // true denotes nonzero protocol. Default for SVG
  const bool useNonzeroRule = true;
  const bool useEvenOddRule = false;
  const bool includeBoundary = false;
  EXPECT_TRUE(in_polygon(PointType {-0.1, 0.0}, poly));
  EXPECT_TRUE(
    in_polygon(PointType {-0.1, 0.0}, poly, includeBoundary, useNonzeroRule));

  // false denotes evenodd protocol
  EXPECT_FALSE(
    in_polygon(PointType {-0.1, 0.0}, poly, includeBoundary, useEvenOddRule));

  // Test in/out on degenerate example in Hormann2001, Figure 6
  axom::Array<PointType> init_vertices(
    {PointType {1, 1}, PointType {1, -2}, PointType {10, -2}, PointType {10, -1}});

  poly = PolygonType(init_vertices);

  axom::Array<PointType> degenerate_vertices({PointType {9, 1},
                                              PointType {8, 0},
                                              PointType {6, 1},
                                              PointType {7, 0},
                                              PointType {5, -1},
                                              PointType {4, 0},
                                              PointType {3, 0}});

  // This point should remain interior when each vertex is added
  EXPECT_TRUE(in_polygon(PointType({2.0, 0.0}), poly));
  for(auto& vertex : degenerate_vertices)
  {
    poly.addVertex(vertex);
    EXPECT_TRUE(in_polygon(PointType({2.0, 0.0}), poly));
  }
}

TEST(primal_polygon, containment_invariants)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  // Invariant to duplicate points
  axom::Array<PointType> vertices(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});
  PolygonType poly;
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 3; j++)  // Duplicate each element 3 times
    {
      poly.addVertex(vertices[i]);
    }
  }

  EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
  EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
  EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));

  // Verify checks up to rotation of vertices
  vertices = axom::Array<PointType>(
    {PointType {0, 0}, PointType {1, 1}, PointType {1, 0}, PointType {0, 1}});

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++)
    {
      poly.addVertex(vertices[(j + i) % 4]);
    }
    EXPECT_TRUE(in_polygon(PointType {0.25, 0.5}, poly));
    EXPECT_TRUE(in_polygon(PointType {0.75, 0.5}, poly));
    EXPECT_FALSE(in_polygon(PointType {0.5, 0.25}, poly));
    EXPECT_FALSE(in_polygon(PointType {0.5, 0.75}, poly));
  }
}

TEST(primal_polygon, containment_edge)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices({PointType {0, 0},
                                   PointType {0.5, 0.5},
                                   PointType {1, 1},
                                   PointType {1, 0},
                                   PointType {0, 1}});
  PolygonType poly(vertices);

  // Edge cases. Default is that points on edges are "outside".
  //  Should work in either in/out protocol
  const bool includeBoundary = true;
  for(bool useNonzeroRule : {true, false})
  {
    EXPECT_TRUE(
      in_polygon(PointType {0, 0.5}, poly, includeBoundary, useNonzeroRule));
    EXPECT_TRUE(
      in_polygon(PointType {0.5, 0.5}, poly, includeBoundary, useNonzeroRule));
    EXPECT_TRUE(
      in_polygon(PointType {.25, .25}, poly, includeBoundary, useNonzeroRule));
    EXPECT_TRUE(
      in_polygon(PointType {.25, .75}, poly, includeBoundary, useNonzeroRule));
    EXPECT_TRUE(
      in_polygon(PointType {1, 0.5}, poly, includeBoundary, useNonzeroRule));
  }

  for(bool useNonzeroRule : {true, false})
  {
    EXPECT_FALSE(
      in_polygon(PointType {0, 0.5}, poly, !includeBoundary, useNonzeroRule));
    EXPECT_FALSE(
      in_polygon(PointType {0.5, 0.5}, poly, !includeBoundary, useNonzeroRule));
    EXPECT_FALSE(
      in_polygon(PointType {.25, .25}, poly, !includeBoundary, useNonzeroRule));
    EXPECT_FALSE(
      in_polygon(PointType {.25, .75}, poly, !includeBoundary, useNonzeroRule));
    EXPECT_FALSE(
      in_polygon(PointType {1, 0.5}, poly, !includeBoundary, useNonzeroRule));
  }

  // Corner cases, where query is on a vertex
  vertices = axom::Array<PointType>({PointType {0, 0},
                                     PointType {0.5, 0},
                                     PointType {1, 0},
                                     PointType {1, 1},
                                     PointType {0, 1}});
  poly = PolygonType(vertices);

  for(auto& vtx : vertices)
  {
    // Nonzero in/out protocol
    EXPECT_TRUE(in_polygon(vtx, poly, includeBoundary, true));
    EXPECT_FALSE(in_polygon(vtx, poly, !includeBoundary, true));

    // Evenodd in/out protocol
    EXPECT_TRUE(in_polygon(vtx, poly, includeBoundary, false));
    EXPECT_FALSE(in_polygon(vtx, poly, !includeBoundary, false));
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, convexity)
{
  using PolygonType = axom::primal::Polygon<double, 2>;
  using PointType = axom::primal::Point<double, 2>;

  axom::Array<PointType> vertices({PointType {0, 0}, PointType {1, 1}});
  PolygonType poly(vertices);

  // Segments and Triangles are always convex
  EXPECT_TRUE(is_convex(poly));

  poly.addVertex(PointType {0, 1});
  EXPECT_TRUE(is_convex(poly));

  axom::Array<PointType> convex_verts = axom::Array<PointType>(
    {PointType {0, 0}, PointType {0, 1}, PointType {1, 1}, PointType {1, 0}});
  axom::Array<PointType> concave_verts = axom::Array<PointType>(
    {PointType {0, 0}, PointType {0, 1}, PointType {0.1, 0.1}, PointType {1, 0}});
  axom::Array<PointType> nonsimple_verts = axom::Array<PointType>(
    {PointType {0, 0}, PointType {1, 1}, PointType {0, 1}, PointType {1, 0}});

  poly.clear();

  // Duplicate points and straight edges should not affect convexity
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 3; j++)  // Duplicate each element 3 times
    {
      poly.addVertex(convex_verts[i]);
    }

    // Add midpoints between each duplicated vertex
    poly.addVertex(
      PointType::midpoint(convex_verts[i], convex_verts[(i + 1) % 4]));
  }

  EXPECT_TRUE(is_convex(poly));

  // Verify checks up to rotation of vertices
  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++)
    {
      poly.addVertex(concave_verts[(j + i) % 4]);
    }
    EXPECT_FALSE(is_convex(poly));
  }

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++)
    {
      poly.addVertex(nonsimple_verts[(j + i) % 4]);
    }
    EXPECT_FALSE(is_convex(poly));
  }

  for(int i = 0; i < 4; i++)
  {
    poly.clear();
    for(int j = 0; j < 4; j++)
    {
      poly.addVertex(convex_verts[(j + i) % 4]);
    }
    EXPECT_TRUE(is_convex(poly));
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, signed_area_2d)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;
  using axom::utilities::abs;

  // test a simple right triangle in CW and CCW orientations
  {
    Polygon2D poly2D_ccw({Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}});
    Polygon2D poly2D_cw({Point2D {0, 0}, Point2D {1, 1}, Point2D {1, 0}});

    // Signed area is positive for CCW and negative for CW
    EXPECT_DOUBLE_EQ(.5, poly2D_ccw.signedArea());
    EXPECT_DOUBLE_EQ(-.5, poly2D_cw.signedArea());

    // The two triangles have reverse orientations (signedArea)
    // but the same (unsigned) area
    EXPECT_DOUBLE_EQ(-poly2D_ccw.signedArea(), poly2D_cw.signedArea());
    EXPECT_DOUBLE_EQ(poly2D_ccw.area(), poly2D_cw.area());

    // compare signed and unsigned areas
    EXPECT_DOUBLE_EQ(poly2D_ccw.area(), poly2D_ccw.signedArea());
    EXPECT_DOUBLE_EQ(-poly2D_cw.area(), poly2D_cw.signedArea());
  }

  // test regular polygons with CW and CCW orienations
  for(int nSides = 3; nSides < 10; ++nSides)
  {
    Polygon2D poly2D_ccw(nSides);
    Polygon2D poly2D_cw(nSides);

    for(int i = 0; i < nSides; ++i)
    {
      const double angle = 2. * M_PI * i / nSides;
      poly2D_ccw.addVertex(Point2D {cos(angle), sin(angle)});
      poly2D_cw.addVertex(Point2D {sin(angle), cos(angle)});
    }

    const double expected_area = nSides / 2. * sin(2 * M_PI / nSides);

    // The areas are the same; signed areas are opposite
    EXPECT_DOUBLE_EQ(expected_area, poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(expected_area, poly2D_cw.area());
    EXPECT_DOUBLE_EQ(poly2D_cw.area(), poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(-poly2D_cw.signedArea(), poly2D_ccw.signedArea());

    EXPECT_DOUBLE_EQ(poly2D_ccw.signedArea(), poly2D_ccw.area());
    EXPECT_DOUBLE_EQ(-poly2D_cw.signedArea(), poly2D_cw.area());
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, area_2d_3d_axis_aligned)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;

  using Polygon3D = axom::primal::Polygon<double, 3>;
  using Point3D = axom::primal::Point<double, 3>;

  // test a simple right triangle
  Polygon2D poly2D({Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}});
  EXPECT_DOUBLE_EQ(.5, poly2D.area());

  // in xy-plane
  Polygon3D poly3Da({Point3D {0, 0, 0}, Point3D {1, 0, 0}, Point3D {1, 1, 0}});
  EXPECT_DOUBLE_EQ(.5, poly3Da.area());

  // in xz-plane
  Polygon3D poly3Db({Point3D {0, 0, 0}, Point3D {1, 0, 0}, Point3D {1, 0, 1}});
  EXPECT_DOUBLE_EQ(.5, poly3Db.area());

  // in yz-plane
  Polygon3D poly3Dc({Point3D {0, 0, 0}, Point3D {0, 1, 0}, Point3D {0, 1, 1}});
  EXPECT_DOUBLE_EQ(.5, poly3Dc.area());
}

//------------------------------------------------------------------------------
TEST(primal_polygon, area_2d_3d_affine_transforms)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;

  using Polygon3D = axom::primal::Polygon<double, 3>;
  using Point3D = axom::primal::Point<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;

  using TransformMatrix = axom::numerics::Matrix<double>;

  // lambda to generate a regular n-sided 2D polygon centered around origin
  auto generateNSidedPolygon = [](int nSides) {
    Polygon2D polygon(nSides);

    for(int i = 0; i < nSides; ++i)
    {
      const double angle = 2. * M_PI * i / nSides;
      polygon.addVertex(Point2D {cos(angle), sin(angle)});
    }

    return polygon;
  };

  // lambda to generate an affine transformation matrix for 2D points
  auto generateTransformMatrix2D =
    [](const Point2D& scale, const Point2D& translate, double rotation_angle) {
      // create scaling matrix
      auto sc_matx = TransformMatrix::identity(3);
      {
        sc_matx(0, 0) = scale[0];
        sc_matx(1, 1) = scale[1];
      }

      // create rotation matrix
      auto rot_matx = TransformMatrix::identity(3);
      {
        const double sinT = std::sin(rotation_angle);
        const double cosT = std::cos(rotation_angle);

        rot_matx(0, 0) = cosT;
        rot_matx(0, 1) = -sinT;
        rot_matx(1, 0) = sinT;
        rot_matx(1, 1) = cosT;
      }

      // create translation matrix
      auto tr_matx = TransformMatrix::identity(3);
      {
        tr_matx(0, 2) = translate[0];
        tr_matx(1, 2) = translate[1];
      }

      // multiply them to get the final transform
      TransformMatrix affine_matx1(3, 3);
      matrix_multiply(rot_matx, sc_matx, affine_matx1);

      TransformMatrix affine_matx2(3, 3);
      matrix_multiply(tr_matx, affine_matx1, affine_matx2);

      EXPECT_NEAR(scale[0] * scale[1], determinant(affine_matx2), EPS);
      return affine_matx2;
    };

  // lambda to generate an affine transformation matrix for 2D points
  auto generateTransformMatrix3D = [](const Point3D& scale,
                                      const Point3D& translate,
                                      const Vector3D& axis,
                                      double angle) {
    // create scaling matrix
    auto sc_matx = TransformMatrix::identity(4);
    {
      sc_matx(0, 0) = scale[0];
      sc_matx(1, 1) = scale[1];
      sc_matx(2, 2) = scale[2];
    }

    // create rotation matrix
    auto rot_matx = TransformMatrix::zeros(4, 4);
    {
      const double sinT = std::sin(angle);
      const double cosT = std::cos(angle);

      const auto unitAxis = axis.unitVector();
      const double& ux = unitAxis[0];
      const double& uy = unitAxis[1];
      const double& uz = unitAxis[2];

      rot_matx(0, 0) = cosT + ux * ux * (1 - cosT);
      rot_matx(0, 1) = ux * uy * (1 - cosT) - uz * sinT;
      rot_matx(0, 2) = ux * uz * (1 - cosT) + uy * sinT;
      rot_matx(1, 0) = uy * ux * (1 - cosT) + uz * sinT;
      rot_matx(1, 1) = cosT + uy * uy * (1 - cosT);
      rot_matx(1, 2) = uy * uz * (1 - cosT) - ux * sinT;
      rot_matx(2, 0) = uz * ux * (1 - cosT) - uy * sinT;
      rot_matx(2, 1) = uz * uy * (1 - cosT) + ux * sinT;
      rot_matx(2, 2) = cosT + uz * uz * (1 - cosT);
      rot_matx(3, 3) = 1;
    }

    // create translation matrix
    auto tr_matx = TransformMatrix::identity(4);
    {
      tr_matx(0, 3) = translate[0];
      tr_matx(1, 3) = translate[1];
      tr_matx(2, 3) = translate[2];
    }

    // multiply them to get the final transform
    TransformMatrix affine_matx1(4, 4);
    matrix_multiply(rot_matx, sc_matx, affine_matx1);
    TransformMatrix affine_matx2(4, 4);
    matrix_multiply(tr_matx, affine_matx1, affine_matx2);

    EXPECT_NEAR(scale[0] * scale[1] * scale[2], determinant(affine_matx2), EPS);
    return affine_matx2;
  };

  // lambda to transform a 2D polygon into 2D
  auto transformedPolygon2d = [](const Polygon2D& poly,
                                 const TransformMatrix& matx) {
    Polygon2D xformed(poly.numVertices());
    for(int i = 0; i < poly.numVertices(); ++i)
    {
      const double vec_in[3] = {poly[i][0], poly[i][1], 1.};
      double vec_out[3] = {0., 0., 0.};
      axom::numerics::matrix_vector_multiply(matx, vec_in, vec_out);
      xformed.addVertex(Point2D {vec_out[0], vec_out[1]});
    }
    return xformed;
  };

  // lambda to transform a 2D polygon into 3D
  auto transformedPolygon3d = [](const Polygon2D& poly,
                                 const TransformMatrix& matx) {
    Polygon3D xformed(poly.numVertices());
    for(int i = 0; i < poly.numVertices(); ++i)
    {
      const double vec_in[4] = {poly[i][0], poly[i][1], 0., 1.};
      double vec_out[4] = {0., 0., 0., 0.};
      axom::numerics::matrix_vector_multiply(matx, vec_in, vec_out);
      xformed.addVertex(Point3D {vec_out[0], vec_out[1], vec_out[2]});
    }
    return xformed;
  };

  const auto scales = axom::Array<double>({-3., -1., -.5, 0., 0.01, 1., 42.3});
  const auto translations = axom::Array<double>({-.5, 0., 1., 42.3});
  const auto angles = axom::Array<double>({-.57, 0., 2. / 3. * M_PI});
  const auto axes = axom::Array<Vector3D>({
    Vector3D {0., 0., 1.},
    Vector3D {0., 1., 0.},
    Vector3D {1., 0., 0.},
    Vector3D {1., 0., 1.},
    Vector3D {1., 1., 1.},
    Vector3D {-2., -5., 0.},
  });

  for(int nSides = 3; nSides < 10; ++nSides)
  {
    Polygon2D polygon2d = generateNSidedPolygon(nSides);
    const double unscaled_area = nSides / 2. * sin(2 * M_PI / nSides);
    EXPECT_DOUBLE_EQ(unscaled_area, polygon2d.area());

    // check area of 2D polygons after affine transforms
    for(double sc_x : scales)
    {
      for(double sc_y : scales)
      {
        for(double tr_x : translations)
        {
          for(double tr_y : translations)
          {
            for(double theta : angles)
            {
              const auto sc = Point2D {sc_x, sc_y};
              const auto tr = Point2D {tr_x, tr_y};
              auto affine_matx = generateTransformMatrix2D(sc, tr, theta);
              auto xformed_polygon = transformedPolygon2d(polygon2d, affine_matx);

              const double expected_area =
                unscaled_area * determinant(affine_matx);
              EXPECT_NEAR(expected_area, xformed_polygon.signedArea(), EPS);

              if(nSides == 3)
              {
                axom::primal::Triangle<double, 2> tri(xformed_polygon[0],
                                                      xformed_polygon[1],
                                                      xformed_polygon[2]);
                EXPECT_NEAR(xformed_polygon.signedArea(), tri.signedArea(), EPS);
              }
            }
          }
        }
      }
    }

    // check area of 3D polygons after affine transforms
    for(double sc_x : scales)
    {
      for(double sc_y : scales)
      {
        for(double tr_x : translations)
        {
          for(double tr_y : translations)
          {
            for(double tr_z : translations)
            {
              for(const auto& axis : axes)
              {
                for(double theta : angles)
                {
                  const auto sc = Point3D {sc_x, sc_y, 1.};
                  const auto tr = Point3D {tr_x, tr_y, tr_z};
                  auto affine_matx =
                    generateTransformMatrix3D(sc, tr, axis, theta);
                  auto xformed_polygon =
                    transformedPolygon3d(polygon2d, affine_matx);

                  const auto expected_area = unscaled_area *
                    axom::utilities::abs(determinant(affine_matx));
                  EXPECT_NEAR(expected_area, xformed_polygon.area(), EPS);

                  if(nSides == 3)
                  {
                    axom::primal::Triangle<double, 3> tri(xformed_polygon[0],
                                                          xformed_polygon[1],
                                                          xformed_polygon[2]);
                    EXPECT_NEAR(xformed_polygon.area(), tri.area(), EPS);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, normal)
{
  using Polygon3D = axom::primal::Polygon<double, 3>;
  using Point3D = axom::primal::Point<double, 3>;

  using Triangle3D = axom::primal::Triangle<double, 3>;
  using Vector3D = axom::primal::Vector<double, 3>;

  // Test a simple right triangle
  {
    Polygon3D poly3D({Point3D {0, 0, 1}, Point3D {1, 0, 0}, Point3D {1, 1, 1}});
    Vector3D poly_normal = poly3D.normal().unitVector();

    Triangle3D tri3D({Point3D {0, 0, 1}, Point3D {1, 0, 0}, Point3D {1, 1, 1}});
    Vector3D tri_normal = tri3D.normal().unitVector();

    for(int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(poly_normal[i], tri_normal[i], 1e-10);
    }
  }

  // Test a planar 5-gon inscribed in a circle
  {
    Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
    Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();
    Vector3D exp_normal = Vector3D::cross_product(v1, v2).unitVector();

    Polygon3D poly(5);

    double angles[5] = {0.0, 0.5, 1.2, 3.0, 5.0};

    for(int i = 0; i < 5; ++i)
    {
      poly.addVertex(Point3D {cos(angles[i]) * v1[0] + sin(angles[i]) * v2[0],
                              cos(angles[i]) * v1[1] + sin(angles[i]) * v2[1],
                              cos(angles[i]) * v1[2] + sin(angles[i]) * v2[2]});
    }

    Vector3D obs_normal = poly.normal().unitVector();

    for(int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(exp_normal[i], obs_normal[i], 1e-10);
    }
  }

  // Test a planar "5-gon" with duplicate vertices and colinear points (robustness)
  {
    Vector3D v1 = Vector3D({0.0, 1.0, 2.0}).unitVector();
    Vector3D v2 = Vector3D({2.0, -1.0, 0.5}).unitVector();
    Vector3D exp_normal = Vector3D::cross_product(v1, v2).unitVector();

    Polygon3D poly(10);

    double angles[8] = {0.0, 0.5, 0.5, 1.2, 3.0, 3.0, 3.0, 5.0};

    // Add point at angle 0
    poly.addVertex(Point3D {cos(angles[0]) * v1[0] + sin(angles[0]) * v2[0],
                            cos(angles[0]) * v1[1] + sin(angles[0]) * v2[1],
                            cos(angles[0]) * v1[2] + sin(angles[0]) * v2[2]});
    // Add a midpoint between angles 0 and 1
    poly.addVertex(Point3D::midpoint(
      poly[0],
      Point3D {cos(angles[1]) * v1[0] + sin(angles[1]) * v2[0],
               cos(angles[1]) * v1[1] + sin(angles[1]) * v2[1],
               cos(angles[1]) * v1[2] + sin(angles[1]) * v2[2]}));

    // Add the rest of the vertices
    for(int i = 1; i < 8; ++i)
    {
      poly.addVertex(Point3D {cos(angles[i]) * v1[0] + sin(angles[i]) * v2[0],
                              cos(angles[i]) * v1[1] + sin(angles[i]) * v2[1],
                              cos(angles[i]) * v1[2] + sin(angles[i]) * v2[2]});
    }

    // Add another midpoint
    poly.addVertex(Point3D::midpoint(poly[0], poly[8]));

    // Verify degeneracies
    EXPECT_NEAR(0.0, Triangle3D(poly[0], poly[1], poly[2]).area(), 1e-10);
    EXPECT_NEAR(0.0, Triangle3D(poly[2], poly[3], poly[4]).area(), 1e-10);

    Vector3D obs_normal = poly.normal().unitVector();

    EXPECT_EQ(10, poly.numVertices());
    for(int i = 0; i < 3; ++i)
    {
      EXPECT_NEAR(exp_normal[i], obs_normal[i], 1e-10);
    }
  }
}

//------------------------------------------------------------------------------
TEST(primal_polygon, reverseOrientation)
{
  using Polygon2D = axom::primal::Polygon<double, 2>;
  using Point2D = axom::primal::Point<double, 2>;

  // Test an odd number of vertices
  {
    Polygon2D poly({Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}});
    poly.reverseOrientation();

    EXPECT_NEAR(poly[0][0], 1, EPS);
    EXPECT_NEAR(poly[0][1], 1, EPS);

    EXPECT_NEAR(poly[1][0], 1, EPS);
    EXPECT_NEAR(poly[1][1], 0, EPS);

    EXPECT_NEAR(poly[2][0], 0, EPS);
    EXPECT_NEAR(poly[2][1], 0, EPS);
  }

  // Test an even number of vertices
  {
    Polygon2D poly(
      {Point2D {0, 0}, Point2D {1, 0}, Point2D {1, 1}, Point2D {0, 1}});
    poly.reverseOrientation();

    EXPECT_NEAR(poly[0][0], 0, EPS);
    EXPECT_NEAR(poly[0][1], 1, EPS);

    EXPECT_NEAR(poly[1][0], 1, EPS);
    EXPECT_NEAR(poly[1][1], 1, EPS);

    EXPECT_NEAR(poly[2][0], 1, EPS);
    EXPECT_NEAR(poly[2][1], 0, EPS);

    EXPECT_NEAR(poly[3][0], 0, EPS);
    EXPECT_NEAR(poly[3][1], 0, EPS);
  }
}

template <typename ExecPolicy>
void check_polygon_policy()
{
  const int NUM_VERTS_SQUARE = 4;

  using Polygon3D =
    axom::primal::Polygon<double, 3, axom::primal::PolygonArray::Static, NUM_VERTS_SQUARE>;
  using Point3D = axom::primal::Point<double, 3>;
  using Polygon2D =
    axom::primal::Polygon<double, 2, axom::primal::PolygonArray::Static, NUM_VERTS_SQUARE>;
  using Point2D = axom::primal::Point<double, 2>;
  using Vector3D = axom::primal::Vector<double, 3>;

  // Get ids of necessary allocators
  const int host_allocator = axom::execution_space<axom::SEQ_EXEC>::allocatorID();
  const int kernel_allocator = axom::execution_space<ExecPolicy>::allocatorID();

  axom::Array<Polygon3D> poly_3d_device(1, 1, kernel_allocator);
  auto poly_3d_view = poly_3d_device.view();
  axom::Array<Polygon2D> poly_2d_device(1, 1, kernel_allocator);
  auto poly_2d_view = poly_2d_device.view();

  axom::Array<Point3D> vertex_mean_3d_device(1, 1, kernel_allocator);
  auto vertex_mean_3d_view = vertex_mean_3d_device.view();
  axom::Array<Point2D> vertex_mean_2d_device(1, 1, kernel_allocator);
  auto vertex_mean_2d_view = vertex_mean_2d_device.view();

  axom::Array<double> area_3d_device(1, 1, kernel_allocator);
  auto area_3d_view = area_3d_device.view();
  axom::Array<double> area_2d_device(1, 1, kernel_allocator);
  auto area_2d_view = area_2d_device.view();

  // 3d only
  axom::Array<Vector3D> normal_3d_device(1, 1, kernel_allocator);
  auto normal_3d_view = normal_3d_device.view();

  axom::for_all<ExecPolicy>(
    1,
    AXOM_LAMBDA(int i) {
      // Initialize to empty polygons
      poly_3d_view[i] = Polygon3D();
      poly_2d_view[i] = Polygon2D();
      poly_3d_view[i].clear();
      poly_2d_view[i].clear();

      // Initialize to triangles
      poly_3d_view[i] = Polygon3D({Point3D({0.0, 0.0, 0.0}),
                                   Point3D({1.0, 0.0, 0.0}),
                                   Point3D({1.0, 1.0, 0.0})});
      poly_2d_view[i] = Polygon2D(
        {Point2D({0.0, 0.0}), Point2D({1.0, 0.0}), Point2D({1.0, 1.0})});

      // Add a vertex to make squares
      (poly_3d_view[i]).addVertex(Point3D({0.0, 1.0, 0.0}));
      (poly_2d_view[i]).addVertex(Point2D({0.0, 1.0}));

      // Collect info about squares
      vertex_mean_3d_view[i] = poly_3d_view[i].vertexMean();
      vertex_mean_2d_view[i] = poly_2d_view[i].vertexMean();
      area_3d_view[i] = poly_3d_view[i].area();
      area_2d_view[i] = poly_2d_view[i].area();
      normal_3d_view[i] = poly_3d_view[i].normal();

      //Sanity check - functions are callable on device
      poly_3d_view[i].numVertices();
      poly_3d_view[i].isValid();
      poly_2d_view[i].numVertices();
      poly_2d_view[i].isValid();

      poly_2d_view[i].reverseOrientation();
      poly_2d_view[i].reverseOrientation();
      poly_3d_view[i].reverseOrientation();
      poly_3d_view[i].reverseOrientation();
    });

  // Copy polygons and data back to host
  axom::Array<Polygon3D> poly_3d_host =
    axom::Array<Polygon3D>(poly_3d_device, host_allocator);
  axom::Array<Polygon2D> poly_2d_host =
    axom::Array<Polygon2D>(poly_2d_device, host_allocator);
  axom::Array<Point3D> vertex_mean_3d_host =
    axom::Array<Point3D>(vertex_mean_3d_device, host_allocator);
  axom::Array<Point2D> vertex_mean_2d_host =
    axom::Array<Point2D>(vertex_mean_2d_device, host_allocator);
  axom::Array<double> area_3d_host =
    axom::Array<double>(area_3d_device, host_allocator);
  axom::Array<double> area_2d_host =
    axom::Array<double>(area_2d_device, host_allocator);
  axom::Array<Vector3D> normal_3d_host =
    axom::Array<Vector3D>(normal_3d_device, host_allocator);

  // Verify values
  EXPECT_EQ(poly_3d_host[0].numVertices(), NUM_VERTS_SQUARE);
  EXPECT_EQ(poly_2d_host[0].numVertices(), NUM_VERTS_SQUARE);

  EXPECT_EQ(vertex_mean_3d_host[0], Point3D({0.5, 0.5, 0}));
  EXPECT_EQ(vertex_mean_2d_host[0], Point2D({0.5, 0.5}));

  EXPECT_DOUBLE_EQ(area_3d_host[0], 1.0);
  EXPECT_DOUBLE_EQ(area_2d_host[0], 1.0);

  EXPECT_EQ(normal_3d_host[0], Vector3D(Point3D({0.0, 0.0, 2.0})));

  EXPECT_TRUE(poly_3d_host[0].isValid());
  EXPECT_TRUE(poly_2d_host[0].isValid());
}

//------------------------------------------------------------------------------
TEST(primal_polygon, polygon_check_seq)
{
  check_polygon_policy<axom::SEQ_EXEC>();
}

#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
  #ifdef AXOM_USE_OPENMP
TEST(primal_polygon, polygon_check_omp)
{
  check_polygon_policy<axom::OMP_EXEC>();
}
  #endif /* AXOM_USE_OPENMP */

  #ifdef AXOM_USE_CUDA
TEST(primal_polygon, polygon_check_cuda)
{
  check_polygon_policy<axom::CUDA_EXEC<256>>();
}
  #endif /* AXOM_USE_CUDA */

  #ifdef AXOM_USE_HIP
TEST(primal_clip, polygon_check_hip)
{
  check_polygon_policy<axom::HIP_EXEC<256>>();
}
  #endif /* AXOM_USE_HIP */

#endif /* AXOM_USE_RAJA && AXOM_USE_UMPIRE */

//------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  result = RUN_ALL_TESTS();

  return result;
}
