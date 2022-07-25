// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file primal_curved_polygon_intersect.cpp
 * \brief This file tests intersections of CurvedPolygon instances
 */

#include "gtest/gtest.h"

#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "axom/fmt.hpp"

#include <numeric>
#include <fstream>

namespace primal = axom::primal;

template <typename CoordType>
void outputAsSVG(
  const std::string& filename,
  const primal::CurvedPolygon<CoordType, 2>& polygon1,
  const primal::CurvedPolygon<CoordType, 2>& polygon2,
  const std::vector<primal::CurvedPolygon<CoordType, 2>> intersectionPolygons)
{
  // Find the bounding box of the set of polygons
  primal::BoundingBox<CoordType, 2> bbox;
  bbox.addBox(polygon1.boundingBox());
  bbox.addBox(polygon2.boundingBox());
  for(const auto& cp : intersectionPolygons)
  {
    bbox.addBox(cp.boundingBox());
  }
  bbox.scale(1.1);

  std::string header = axom::fmt::format(
    "<svg viewBox='{} {} {} {}' xmlns='http://www.w3.org/2000/svg'>",
    bbox.getMin()[0],
    bbox.getMin()[1],
    bbox.range()[0],
    bbox.range()[1]);
  std::string footer = "</svg>";

  // lambda to convert a CurvedPolygon to an SVG path string
  auto cpToSVG = [](const primal::CurvedPolygon<CoordType, 2>& cp) {
    axom::fmt::memory_buffer out;
    bool is_first = true;

    for(auto& curve : cp.getEdges())
    {
      // Only write out first point for first edge
      if(is_first)
      {
        axom::fmt::format_to(out, "M {} {} ", curve[0][0], curve[0][1]);
        is_first = false;
      }

      switch(curve.getOrder())
      {
      case 1:
        axom::fmt::format_to(out, "L {} {} ", curve[1][0], curve[1][1]);
        break;
      case 2:
        axom::fmt::format_to(out,
                             "Q {} {}, {} {} ",
                             curve[1][0],
                             curve[1][1],
                             curve[2][0],
                             curve[2][1]);
        break;
      case 3:
        axom::fmt::format_to(out,
                             "C {} {}, {} {}, {} {} ",
                             curve[1][0],
                             curve[1][1],
                             curve[2][0],
                             curve[2][1],
                             curve[3][0],
                             curve[3][1]);
        break;
      default:
        SLIC_WARNING(
          "Unsupported case: can only output up to cubic curves as SVG.");
      }
    }
    return axom::fmt::format("    <path d='{} Z' />\n",
                             axom::fmt::to_string(out));
  };

  std::string poly1Group;
  std::string poly2Group;
  std::string intersectionGroup;

  // render polygon1 as SVG
  {
    axom::fmt::memory_buffer out;
    axom::fmt::format_to(
      out,
      "  <g id='source_mesh' stroke='black' stroke-width='.01' "
      "fill='red' fill-opacity='.7'>\n");
    axom::fmt::format_to(out, cpToSVG(polygon1));
    axom::fmt::format_to(out, "  </g>\n");
    poly1Group = axom::fmt::to_string(out);
  }

  // render polygon2 as SVG
  {
    axom::fmt::memory_buffer out;
    axom::fmt::format_to(
      out,
      "  <g id='target_mesh' stroke='black' stroke-width='.01' "
      "fill='blue' fill-opacity='.7'>\n");
    axom::fmt::format_to(out, cpToSVG(polygon2));
    axom::fmt::format_to(out, "  </g>\n");
    poly2Group = axom::fmt::to_string(out);
  }

  //render intersection polygons as SVG
  {
    axom::fmt::memory_buffer out;
    axom::fmt::format_to(
      out,
      "  <g id='intersection_mesh' stroke='black' stroke-width='.01' "
      "fill='green' fill-opacity='.7'>\n");
    for(auto& cp : intersectionPolygons)
    {
      axom::fmt::format_to(out, cpToSVG(cp));
    }
    axom::fmt::format_to(out, "  </g>\n");
    intersectionGroup = axom::fmt::to_string(out);
  }

  // Write the file
  {
    std::ofstream fs(filename);
    fs << header << std::endl;
    fs << poly1Group << std::endl;
    fs << poly2Group << std::endl;
    fs << intersectionGroup << std::endl;
    fs << footer << std::endl;
  }
}

/*!
 * Helper function to compute the set of intersection polygons given two input polygons
 * and to check that they match expectations, stored in \a expbPolygon.
 * Intersection polygon is computed to within tolerance \a eps and checks use \a test_eps.
 */
template <typename CoordType>
void checkIntersection(
  const primal::CurvedPolygon<CoordType, 2>& bPolygon1,
  const primal::CurvedPolygon<CoordType, 2>& bPolygon2,
  const std::vector<primal::CurvedPolygon<CoordType, 2>> expbPolygon,
  const double eps = 1e-15,
  const double test_eps = 1e-13)
{
  constexpr int DIM = 2;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using BezierCurveType = typename CurvedPolygonType::BezierCurveType;
  using PointType = typename BezierCurveType::PointType;

  std::vector<CurvedPolygonType> intersectionPolys;

  //Compute intersection using algorithm with tolerance of eps
  intersect(bPolygon1, bPolygon2, intersectionPolys, eps);
  //Check that expected number of intersection regions are found
  ASSERT_EQ(expbPolygon.size(), intersectionPolys.size());

  //Check that expected intersection curves are found to within test_eps
  const int nPolygons = expbPolygon.size();
  for(int p = 0; p < nPolygons; ++p)
  {
    const CurvedPolygonType& polyExp = expbPolygon[p];
    const CurvedPolygonType& polyActual = intersectionPolys[p];
    EXPECT_EQ(polyExp.numEdges(), polyActual.numEdges());

    const int nEdges = polyExp.numEdges();
    for(int e = 0; e < nEdges; ++e)
    {
      const BezierCurveType& curveExp = polyExp[e];
      const BezierCurveType& curveActual = polyActual[e];
      EXPECT_EQ(curveExp.getOrder(), curveActual.getOrder());

      const int nPts = curveExp.getOrder() + 1;
      for(int idx = 0; idx < nPts; ++idx)
      {
        const PointType& ptExp = curveExp[idx];
        const PointType& ptActual = curveActual[idx];

        for(int d = 0; d < DIM; ++d)
        {
          EXPECT_NEAR(ptExp[d], ptActual[d], test_eps)
            << "Difference in polygon " << p << " edge " << e
            << " control point " << idx << " dimension " << d;
        }
      }
    }
  }
}

/*!
 * Helper function to create a CurvedPolygon from a list of control points and a list
 * of orders of component curves. Control points should be given as a list of Points
 * in order of orientation with no duplicates except that the first control point
 * should also be the last control point (if the polygon is closed).  Orders should
 * be given as a list of ints in order of orientation, representing the orders of the component curves.
 */
template <typename CoordType>
primal::CurvedPolygon<CoordType, 2> createPolygon(
  const std::vector<primal::Point<CoordType, 2>> ControlPoints,
  const std::vector<int> orders)
{
  using PointType = primal::Point<CoordType, 2>;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, 2>;
  using BezierCurveType = primal::BezierCurve<CoordType, 2>;

  const int num_edges = orders.size();
  const int num_unique_control_points = ControlPoints.size();

  //checks if the orders and control points given will yield a valid curved polygon
  {
    const int sum_of_orders = std::accumulate(orders.begin(), orders.end(), 0);
    EXPECT_EQ(sum_of_orders + 1, num_unique_control_points);
  }

  //Converts the control points to BezierCurves of specified orders and stores them in a CurvedPolygon object.
  CurvedPolygonType bPolygon;
  int iter = 0;
  for(int j = 0; j < num_edges; ++j)
  {
    std::vector<PointType> subCP;
    subCP.assign(ControlPoints.begin() + iter,
                 ControlPoints.begin() + iter + orders[j] + 1);
    BezierCurveType addCurve(subCP, orders[j]);
    bPolygon.addEdge(addCurve);
    iter += (orders[j]);
  }
  return bPolygon;
}

//----------------------------------------------------------------------------------

TEST(primal_curvedpolygon, detail_intersection_type)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  using primal::detail::getJunctionIntersectionType;
  using primal::detail::JunctionIntersectionType;

  SLIC_INFO("Tests various junction intersections");

  // the first pair is a horizontal line from left to right
  std::vector<int> orders = {1, 1};
  std::vector<PointType> CP = {PointType {-1, 0},
                               PointType {0, 0},
                               PointType {1, 0}};
  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  // A type2 case
  {
    std::vector<PointType> CP2 = {PointType {1, -1},
                                  PointType {0, 0},
                                  PointType {-1, -1}};
    CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

    auto xType = getJunctionIntersectionType(bPolygon1[0],
                                             bPolygon1[1],
                                             bPolygon2[0],
                                             bPolygon2[1]);
    EXPECT_EQ(JunctionIntersectionType::Type2, xType);
  }
}

//----------------------------------------------------------------------------------

TEST(primal_curvedpolygon, intersection_squares)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two squares (single region)");

  std::vector<int> orders = {1, 1, 1, 1};

  // Unit square scaled by 2
  std::vector<PointType> CP = {PointType {0, 0},
                               PointType {2, 0},
                               PointType {2, 2},
                               PointType {0, 2},
                               PointType {0, 0}};

  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  // Unit square scaled by 2 and offset by (-1,-1)
  std::vector<PointType> CP2 = {PointType {-1, -1},
                                PointType {1, -1},
                                PointType {1, 1},
                                PointType {-1, 1},
                                PointType {-1, -1}};

  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

  // Intersection should be a unit square
  std::vector<PointType> expCP = {PointType {1, 0},
                                  PointType {1, 1},
                                  PointType {0, 1},
                                  PointType {0, 0},
                                  PointType {1, 0}};
  std::vector<int> exporders = {1, 1, 1, 1};
  CurvedPolygonType expbPolygon = createPolygon(expCP, exporders);

  std::vector<CurvedPolygonType> expbPolygons {expbPolygon};
  checkIntersection(bPolygon1, bPolygon2, expbPolygons);

  // Output intersections as SVG
  {
    std::vector<CurvedPolygonType> intersections;
    intersect(bPolygon1, bPolygon2, intersections, 1e-15);
    outputAsSVG("curved_polygon_intersections_linear_squares.svg",
                bPolygon1,
                bPolygon2,
                intersections);
  }
}

//----------------------------------------------------------------------------------

TEST(primal_curvedpolygon, intersection_triangle_linear)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two linear triangles (single region)");

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.3, 2.0},
                               PointType {0.0, 1.6},
                               PointType {0.6, 1.2}};
  std::vector<int> orders = {1, 1, 1};

  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  std::vector<PointType> CP2 = {PointType {0.71, 1.31},
                                PointType {0.41, 2.11},
                                PointType {0.11, 1.71},
                                PointType {0.71, 1.31}};

  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

  std::vector<PointType> expCP = {
    PointType {0.3091666666666666666666, 1.9755555555555555555},
    PointType {0.11, 1.71},
    PointType {0.5083333333333333333, 1.44444444444444444444},
    PointType {0.3091666666666666666666, 1.9755555555555555555}};
  std::vector<int> exporders = {1, 1, 1};
  CurvedPolygonType expbPolygon = createPolygon(expCP, exporders);

  std::vector<CurvedPolygonType> expbPolygons {expbPolygon};
  checkIntersection(bPolygon1, bPolygon2, expbPolygons);

  // Output intersections as SVG
  {
    std::vector<CurvedPolygonType> intersections;
    intersect(bPolygon1, bPolygon2, intersections, 1e-15);
    outputAsSVG("curved_polygon_intersections_linear_triangles.svg",
                bPolygon1,
                bPolygon2,
                intersections);
  }
}

//----------------------------------------------------------------------------------

TEST(primal_curvedpolygon, intersections_triangle_rectangle)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO(
    "Test several intersection cases b/w a linear triangle and rectangle");

  // Rectangle with bounds, -5 <= x < 5 ; and 0 <= y <= 2
  std::vector<PointType> rectanglePts = {PointType {-5, 0},
                                         PointType {5, 0},
                                         PointType {5, 2},
                                         PointType {-5, 2},
                                         PointType {-5, 0}};
  std::vector<int> rectangleOrders = {1, 1, 1, 1};

  // Equilateral triangle with base from -3 <= x <= 3 and height 1
  CurvedPolygonType bRectangle = createPolygon(rectanglePts, rectangleOrders);

  std::vector<PointType> triPts = {PointType {-3, -2},
                                   PointType {3, -2},
                                   PointType {0, -1},
                                   PointType {-3, -2}};
  std::vector<int> triOrders = {1, 1, 1};

  const bool bVerbose = true;
  const double EPS = 1e-7;
  const double SQ_EPS = EPS * EPS;

  // No intersection case: Triangle is below the rectangle
  {
    CurvedPolygonType bTriangle = createPolygon(triPts, triOrders);
    primal::detail::DirectionalWalk<double> walk(bVerbose);
    int nIntersections =
      walk.splitPolygonsAlongIntersections(bRectangle, bTriangle, SQ_EPS);
    EXPECT_EQ(0, nIntersections);

    std::vector<CurvedPolygonType> intersections;
    EXPECT_FALSE(primal::intersect(bRectangle, bTriangle, intersections, EPS));

    outputAsSVG("intersection_test_tri_rect_none.svg",
                bRectangle,
                bTriangle,
                intersections);
  }

  // No intersection case: Triangle apex grazes rectangle
  {
    triPts[2][1] = 0;
    CurvedPolygonType bTriangle = createPolygon(triPts, triOrders);

    primal::detail::DirectionalWalk<double> walk(bVerbose);
    int nIntersections =
      walk.splitPolygonsAlongIntersections(bRectangle, bTriangle, SQ_EPS);
    EXPECT_EQ(1, nIntersections);

    std::vector<CurvedPolygonType> walkIntersections;
    walk.findIntersectionRegions(walkIntersections);

    // EXPECT_EQ(0, walkIntersections.size());   // Warning: Not properly handled yet

    outputAsSVG("intersection_test_tri_rect_lower_graze.svg",
                bRectangle,
                bTriangle,
                walkIntersections);

    std::vector<CurvedPolygonType> intersections;
    // Warning: Not properly handled yet
    // EXPECT_FALSE(
    primal::intersect(bRectangle, bTriangle, intersections, EPS)
      //)
      ;
  }

  // Simple intersection case: Triangle apex inside rectangle
  {
    triPts[2][1] = 1;
    CurvedPolygonType bTriangle = createPolygon(triPts, triOrders);

    primal::detail::DirectionalWalk<double> walk(bVerbose);
    int nIntersections =
      walk.splitPolygonsAlongIntersections(bRectangle, bTriangle, SQ_EPS);
    EXPECT_EQ(2, nIntersections);

    std::vector<CurvedPolygonType> walkIntersections;
    walk.findIntersectionRegions(walkIntersections);
    EXPECT_EQ(1, walkIntersections.size());

    outputAsSVG("intersection_test_tri_rect_intersect_tri.svg",
                bRectangle,
                bTriangle,
                walkIntersections);

    std::vector<PointType> expCP = {PointType {1, 0},
                                    PointType {0, 1},
                                    PointType {-1, 0},
                                    PointType {1, 0}};
    std::vector<int> exporders = {1, 1, 1};

    std::vector<CurvedPolygonType> expbPolygons = {
      createPolygon(expCP, exporders)};

    checkIntersection(bRectangle, bTriangle, expbPolygons);
  }

  // Grazing intersection case: Triangle apex intersects top
  {
    triPts[2][1] = 2;
    CurvedPolygonType bTriangle = createPolygon(triPts, triOrders);

    primal::detail::DirectionalWalk<double> walk(bVerbose);
    int nIntersections =
      walk.splitPolygonsAlongIntersections(bRectangle, bTriangle, SQ_EPS);
    EXPECT_EQ(3, nIntersections);

    std::vector<CurvedPolygonType> walkIntersections;
    walk.findIntersectionRegions(walkIntersections);
    EXPECT_EQ(1, walkIntersections.size());

    //EXPECT_EQ(3, walkIntersections[0].numEdges()); // Warning: Not properly handled yet

    outputAsSVG("intersection_test_tri_rect_upper_graze.svg",
                bRectangle,
                bTriangle,
                walkIntersections);

    std::vector<CurvedPolygonType> intersections;
    // Warning: Not properly handled yet
    EXPECT_TRUE(primal::intersect(bRectangle, bTriangle, intersections, EPS));
  }

  // intersection case: Triangle apex above rectangle
  {
    triPts[2][1] = 3;
    CurvedPolygonType bTriangle = createPolygon(triPts, triOrders);

    primal::detail::DirectionalWalk<double> walk(bVerbose);
    int nIntersections =
      walk.splitPolygonsAlongIntersections(bRectangle, bTriangle, SQ_EPS);
    EXPECT_EQ(4, nIntersections);

    std::vector<CurvedPolygonType> walkIntersections;
    walk.findIntersectionRegions(walkIntersections);
    EXPECT_EQ(1, walkIntersections.size());

    outputAsSVG("intersection_test_tri_rect_intersect_rect.svg",
                bRectangle,
                bTriangle,
                walkIntersections);

    std::vector<PointType> expCP = {PointType {1.8, 0},
                                    PointType {0.6, 2},
                                    PointType {-0.6, 2},
                                    PointType {-1.8, 0},
                                    PointType {1.8, 0}};
    std::vector<int> exporders = {1, 1, 1, 1};

    std::vector<CurvedPolygonType> expbPolygons = {
      createPolygon(expCP, exporders)};

    checkIntersection(bRectangle, bTriangle, expbPolygons);
  }
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersecting two quadratic triangles (single region)");

  CurvedPolygonType bPolygon;
  EXPECT_EQ(0, bPolygon.numEdges());

  std::vector<PointType> CP = {PointType {0.6, 1.2},
                               PointType {0.4, 1.3},
                               PointType {0.3, 2.0},
                               PointType {0.27, 1.5},
                               PointType {0.0, 1.6},
                               PointType {0.1, 1.5},
                               PointType {0.6, 1.2}};
  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP, orders);

  std::vector<PointType> CP2 = {PointType {0.71, 1.31},
                                PointType {0.51, 1.41},
                                PointType {0.41, 2.11},
                                PointType {0.38, 1.61},
                                PointType {0.11, 1.71},
                                PointType {0.21, 1.61},
                                PointType {0.71, 1.31}};
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

  std::vector<PointType> expCP = {
    PointType {0.335956890729522, 1.784126953773395},
    PointType {0.297344765794753, 1.718171485335525},
    PointType {0.2395677533016981, 1.700128235793371},
    PointType {0.221884203146682, 1.662410644580941},
    PointType {0.199328465398189, 1.636873522352205},
    PointType {0.277429214338182, 1.579562422716502},
    PointType {0.408882616650578, 1.495574996394597},
    PointType {0.368520120719339, 1.616453177259694},
    PointType {0.335956890729522, 1.784126953773394}};
  std::vector<int> exporders = {2, 2, 2, 2};
  CurvedPolygonType expbPolygon = createPolygon(expCP, exporders);
  std::vector<CurvedPolygonType> expbPolygons = {expbPolygon};

  checkIntersection(bPolygon1, bPolygon2, expbPolygons, 1e-15, 1e-13);
}

//----------------------------------------------------------------------------------
TEST(primal_curvedpolygon, intersection_triangle_quadratic_two_regions)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two quadratic triangles (two regions)");

  std::vector<PointType> CP1 = {PointType {0.6, 1.2},
                                PointType {0.4, 1.3},
                                PointType {0.3, 2.0},
                                PointType {0.27, 1.5},
                                PointType {0.0, 1.6},
                                PointType {0.1, 1.5},
                                PointType {0.6, 1.2}};

  std::vector<PointType> CP2 = {PointType {1.0205, 1.6699},
                                PointType {0.8339, 1.5467},
                                PointType {0.1777, 1.8101},
                                PointType {0.5957, 1.5341},
                                PointType {0.3741, 1.3503},
                                PointType {0.5107, 1.3869},
                                PointType {1.0205, 1.6699}};
  std::vector<int> orders = {2, 2, 2};

  std::vector<PointType> expCP1 = {
    PointType {0.343364196589264, 1.747080669655736},
    PointType {0.305984025190458, 1.760433098612141},
    PointType {0.266743999290327, 1.775316659915674},
    PointType {0.263419346128088, 1.763343410502168},
    PointType {0.259796003065908, 1.752116885838515},
    PointType {0.320641367919239, 1.705796408318085},
    PointType {0.362111919147859, 1.662268860466508},
    PointType {0.352450139541348, 1.702947255097842},
    PointType {0.343364196589264, 1.747080669655736}};

  std::vector<PointType> expCP2 = {
    PointType {0.435276730907216, 1.423589798138227},
    PointType {0.416268597450954, 1.385275578571685},
    PointType {0.374100000000000, 1.350300000000000},
    PointType {0.404839872482010, 1.358536305511285},
    PointType {0.454478985809487, 1.379250566393211},
    PointType {0.444689566319939, 1.400290430035245},
    PointType {0.435276730907216, 1.423589798138227}};

  std::vector<int> exporder1 = {2, 2, 2, 2};
  std::vector<int> exporder2 = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);
  CurvedPolygonType expbPolygon1 = createPolygon(expCP1, exporder1);
  CurvedPolygonType expbPolygon2 = createPolygon(expCP2, exporder2);
  std::vector<CurvedPolygonType> expIntersections = {expbPolygon1, expbPolygon2};

  checkIntersection(bPolygon1, bPolygon2, expIntersections);

  // Output intersections as SVG
  {
    std::vector<CurvedPolygonType> intersections;
    intersect(bPolygon1, bPolygon2, intersections, 1e-15);
    outputAsSVG(
      "curved_polygon_intersections_quadratic_triangles_two_regions.svg",
      bPolygon1,
      bPolygon2,
      intersections);
  }
}

TEST(primal_curvedpolygon, area_intersection_triangle_inclusion)
{
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of two quadratic triangles (inclusion)");

  std::vector<PointType> CP1 = {PointType {0.0, 0.0},
                                PointType {0.5, 0.0},
                                PointType {1.0, 0.0},
                                PointType {0.5, 0.5},
                                PointType {0.0, 1.0},
                                PointType {0.0, 0.5},
                                PointType {0.0, 0.0}};

  std::vector<PointType> CP2 = {PointType {0.05, 0.05},
                                PointType {0.30, 0.05},
                                PointType {0.55, 0.05},
                                PointType {0.30, 0.30},
                                PointType {0.05, 0.55},
                                PointType {0.05, 0.30},
                                PointType {0.05, 0.05}};
  std::vector<int> orders = {2, 2, 2};
  CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
  CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);
  std::vector<CurvedPolygonType> expIntersections = {bPolygon2};

  checkIntersection(bPolygon1, bPolygon2, expIntersections);
  checkIntersection(bPolygon2, bPolygon1, expIntersections);
}

TEST(primal_curvedpolygon, doubleIntersection)
{
  const double EPS = 1e-8;
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Tests multiple intersections along an edge");

  {
    // Unit square
    std::vector<PointType> CP1 = {PointType {0, 0},
                                  PointType {1, 0},
                                  PointType {1, 1},
                                  PointType {0, 1},
                                  PointType {0, 0}};
    std::vector<int> orders1 = {1, 1, 1, 1};

    // Bi-gon defined by a quadratic edge and a straight line
    std::vector<PointType> CP2 = {PointType {0.8, .25},
                                  PointType {2.0, .50},
                                  PointType {0.8, .75},
                                  PointType {0.8, .25}};
    std::vector<int> orders2 = {2, 1};

    CurvedPolygonType bPolygon1 = createPolygon(CP1, orders1);
    CurvedPolygonType bPolygon2 = createPolygon(CP2, orders2);

    // Find intersections using DirectionalWalk class
    double walkIntersectionArea = 0.;
    {
      const bool verbose = true;
      primal::detail::DirectionalWalk<double> walk(verbose);
      int numIntersections =
        walk.splitPolygonsAlongIntersections(bPolygon1, bPolygon2, EPS * EPS);

      EXPECT_EQ(2, numIntersections);
      EXPECT_NEAR(primal::area(bPolygon1), primal::area(walk.psplit[0]), EPS);
      EXPECT_NEAR(primal::area(bPolygon2), primal::area(walk.psplit[1]), EPS);

      SLIC_INFO("Checking for intersections between " << bPolygon1 << " and "
                                                      << bPolygon2);

      std::vector<CurvedPolygonType> regions;
      walk.findIntersectionRegions(regions);
      EXPECT_EQ(1, regions.size());

      if(!regions.empty())
      {
        EXPECT_EQ(4, regions[0].numEdges());
        walkIntersectionArea = primal::area(regions[0]);
      }
      else
      {
        FAIL() << "Expected the two polygons to intersect";
      }

      SLIC_INFO("Found intersections (directional walk): ");
      for(auto cp : regions)
      {
        SLIC_INFO("\t" << cp);
      }
    }

    // Find interesections using primal::intersect
    double directIntersectionArea = 0.;
    {
      std::vector<CurvedPolygonType> regions;
      bool intersects = intersect(bPolygon1, bPolygon2, regions, EPS);
      EXPECT_TRUE(intersects);
      EXPECT_EQ(1, regions.size());

      if(!regions.empty())
      {
        EXPECT_EQ(4, regions[0].numEdges());
        directIntersectionArea = primal::area(regions[0]);
      }
      else
      {
        FAIL() << "Expected the two polygons to intersect";
      }

      SLIC_INFO("Found intersections (direct): ");
      for(auto cp : regions)
      {
        SLIC_INFO("\t" << cp);
      }
    }

    EXPECT_NEAR(walkIntersectionArea, directIntersectionArea, EPS);
  }
}

TEST(primal_curvedpolygon, regression)
{
  const double EPS = 1e-8;
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of pairs of polygons from regression data");

  // First test: Intersecting a pair of linear quadrilaterals
  // Note: One of the intersections happens betwen a vertex and edge
  // and is not yet properly handled by primal::intersect
  {
    std::vector<PointType> CP1 = {PointType {1, -2},
                                  PointType {0, -1},
                                  PointType {0, 1},
                                  PointType {1, 2},
                                  PointType {1, -2}};

    std::vector<PointType> CP2 = {PointType {-.9, -2},
                                  PointType {0.1, -1},
                                  PointType {2.1, -1},
                                  PointType {3.1, -2},
                                  PointType {-.9, -2}};
    std::vector<int> orders = {1, 1, 1, 1};

    CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
    CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);
    std::vector<CurvedPolygonType> expIntersections;

    intersect(bPolygon1, bPolygon2, expIntersections, EPS);

    SLIC_INFO("There were " << expIntersections.size()
                            << " intersection polygons");
    for(auto cp : expIntersections)
    {
      SLIC_INFO("\t" << cp);
    }
  }
}

TEST(primal_curvedpolygon, adjacent_squares)
{
  const double EPS = 1e-8;
  const int DIM = 2;
  using CoordType = double;
  using CurvedPolygonType = primal::CurvedPolygon<CoordType, DIM>;
  using PointType = primal::Point<CoordType, DIM>;

  SLIC_INFO("Test intersection of pairs of adjacent squares");

  // Sharing an edge
  {
    std::vector<PointType> CP1 = {PointType {-1, 0},
                                  PointType {0, 0},
                                  PointType {0, 1},
                                  PointType {-1, 1},
                                  PointType {-1, 0}};

    std::vector<PointType> CP2 = {PointType {0, 0},
                                  PointType {0, 1},
                                  PointType {1, 1},
                                  PointType {1, 0},
                                  PointType {0, 0}};
    std::vector<int> orders = {1, 1, 1, 1};

    CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
    CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

    {
      std::vector<CurvedPolygonType> expIntersections;
      intersect(bPolygon1, bPolygon2, expIntersections, EPS);

      EXPECT_TRUE(expIntersections.empty());

      if(!expIntersections.empty())
      {
        SLIC_INFO("There were " << expIntersections.size()
                                << " intersection polygons");
        for(auto cp : expIntersections)
        {
          SLIC_INFO("\t" << cp);
        }
        outputAsSVG("adjacent_squares_edge_LR.svg",
                    bPolygon1,
                    bPolygon2,
                    expIntersections);
      }
    }

    {
      std::vector<CurvedPolygonType> expIntersections;
      intersect(bPolygon2, bPolygon1, expIntersections, EPS);

      EXPECT_TRUE(expIntersections.empty());

      if(!expIntersections.empty())
      {
        SLIC_INFO("There were " << expIntersections.size()
                                << " intersection polygons");
        for(auto cp : expIntersections)
        {
          SLIC_INFO("\t" << cp);
        }
        outputAsSVG("adjacent_squares_edge_RL.svg",
                    bPolygon2,
                    bPolygon1,
                    expIntersections);
      }
    }
  }

  // Sharing a vertex
  if(false)
  {
    std::vector<PointType> CP1 = {PointType {-1, 0},
                                  PointType {0, 0},
                                  PointType {0, 1},
                                  PointType {-1, 1},
                                  PointType {-1, 0}};

    std::vector<PointType> CP2 = {PointType {0, 1},
                                  PointType {0, 2},
                                  PointType {1, 2},
                                  PointType {1, 1},
                                  PointType {0, 1}};
    std::vector<int> orders = {1, 1, 1, 1};

    CurvedPolygonType bPolygon1 = createPolygon(CP1, orders);
    CurvedPolygonType bPolygon2 = createPolygon(CP2, orders);

    {
      std::vector<CurvedPolygonType> expIntersections;
      intersect(bPolygon1, bPolygon2, expIntersections, EPS);

      EXPECT_TRUE(expIntersections.empty());

      if(!expIntersections.empty())
      {
        SLIC_INFO("There were " << expIntersections.size()
                                << " intersection polygons");
        for(auto cp : expIntersections)
        {
          SLIC_INFO("\t" << cp);
        }
        outputAsSVG("adjacent_squares_corner_LR.svg",
                    bPolygon1,
                    bPolygon2,
                    expIntersections);
      }
    }

    if(false)
    {
      std::vector<CurvedPolygonType> expIntersections;
      intersect(bPolygon2, bPolygon1, expIntersections, EPS);

      EXPECT_TRUE(expIntersections.empty());

      if(!expIntersections.empty())
      {
        SLIC_INFO("There were " << expIntersections.size()
                                << " intersection polygons");
        for(auto cp : expIntersections)
        {
          SLIC_INFO("\t" << cp);
        }
        outputAsSVG("adjacent_squares_corner_RL.svg",
                    bPolygon2,
                    bPolygon1,
                    expIntersections);
      }
    }
  }
}

//----------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  axom::slic::SimpleLogger logger;

  result = RUN_ALL_TESTS();

  return result;
}
