// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_DISCRETIZE_TEST_HPP_
#define QUEST_DISCRETIZE_TEST_HPP_

#include "axom/quest/Discretize.hpp"  // quest::Discretize
#include "axom/config.hpp"
#include "axom/slic/interface/slic.hpp"
#include "axom/core.hpp"
#include "axom/primal.hpp"

using SphereType = axom::primal::Sphere<double, 3>;
using OctType = axom::primal::Octahedron<double, 3>;
using Point3D = axom::primal::Point<double, 3>;
using Point2D = axom::primal::Point<double, 2>;
using NAType = axom::primal::NumericArray<double, 3>;

// gtest includes
#include "gtest/gtest.h"


#include <vector>

//------------------------------------------------------------------------------
bool check_generation(OctType *& standard,
                      OctType *& test,
                      int generation,
                      int offset,
                      int count)
{
  std::vector<int> matched(count, 0);

  for(int i = 0; i < count; ++i)
  {
    int has_matched = -1;
    for(int j = 0; j < count; ++j)
    {
      if(!matched[j] && standard[offset + i].equals(test[offset + j]))
      {
        matched[j] = 1;
        has_matched = offset + j;
      }
    }
    if(has_matched < 0)
    {
      std::cout << "Gen " << generation << " standard oct " << (offset + i)
                << " didn't match: " << standard[offset + i] << std::endl;
    }
  }

  int matchcount = 0;
  for(int i = 0; i < count; ++i)
  {
    matchcount += matched[i];
  }

  bool matches = (matchcount == count);

  if(!matches)
  {
    // a little more reporting here
    std::cout << "Generation " << generation << " had " << (count - matchcount)
              << " test octahedra not matched to standard octahedra:"
              << std::endl;
    for(int i = 0; i < count; ++i)
    {
      if(!matched[i])
      {
        std::cout << "Test oct " << (offset + i)
                  << " not matched:" << test[offset + i] << std::endl;
      }
    }
  }

  return matches;
}


/* Return a handwritten list of the octahedra discretizing a one-segment polyline.
 *
 * The routine allocates and returns an array of octahedrons pointed to by out.
 * The caller must free that array.
 */
void discretized_segment(Point2D a, Point2D b, OctType *& out)
{
  // We're going to return three generations in the out-vector:
  // one in the first generation, three in the second (covering each of the
  // side-walls of the first generation), eight in the third (covering the
  // two exposed side-walls in each of the second-gen octs).
  //
  // These octahedra are actually triangular prisms (if the polyline
  // is parallel to the X-axis) or truncated tetrahedra.  They have two parallel
  // triangle end-caps and three quadrilateral side-walls.  Because we store
  // them in Octahedron objects, those side-walls are actually two triangles
  // that are coplanar.
  constexpr int FIRST_GEN_COUNT = 1;
  constexpr int SECOND_GEN_COUNT = 3;
  constexpr int THIRD_GEN_COUNT = 6;
  out = axom::allocate<OctType>(FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT);

  // The first generation puts a triangle in the end-discs of the truncated
  // cones with vertices at 12, 4, and 8 o'clock.
  // The second generation puts three triangles in the end-discs;
  //     the first has vertices at 12, 2, and 4 o'clock,
  //     the second at 4, 6, and 8 o'clock,
  //     the third at 8, 10, and 12 o'clock.
  // The third generation puts six more triangles in the ends;
  //     first at 12, 1, 2 o'clock,
  //     second at 2, 3, 4 o'clock,
  //     third at 4, 5, 6, o'clock,
  //     fourth at 6, 7, 8 o'clock,
  //     fifth at 8, 9, 10 o'clock,
  //     sixth at 10, 11, 12 o'clock.
  // These are constructed with 1, 1/2, and sqrt(3)/2 in discs perpendicular
  // to the X-axis at a[0] and b[0].  We're looking out the X-axis so a
  // clockwise spin corresponds to the right-hand spin in the function
  // from_segment() in Discretize.cpp.

  const double SQ_3_2 = sqrt(3.) / 2.;

  // clang-format off
  Point3D Ia   {a[0],    -0.5*a[1],  SQ_3_2*a[1]};
  Point3D IIa  {a[0], -SQ_3_2*a[1],     0.5*a[1]};
  Point3D IIIa {a[0],      -1*a[1],       0*a[1]};
  Point3D IVa  {a[0], -SQ_3_2*a[1],    -0.5*a[1]};
  Point3D Va   {a[0],    -0.5*a[1], -SQ_3_2*a[1]};
  Point3D VIa  {a[0],       0*a[1],      -1*a[1]};
  Point3D VIIa {a[0],     0.5*a[1], -SQ_3_2*a[1]};
  Point3D VIIIa{a[0],  SQ_3_2*a[1],    -0.5*a[1]};
  Point3D IXa  {a[0],       1*a[1],       0*a[1]};
  Point3D Xa   {a[0],  SQ_3_2*a[1],     0.5*a[1]};
  Point3D XIa  {a[0],     0.5*a[1],  SQ_3_2*a[1]};
  Point3D XIIa {a[0],       0*a[1],       1*a[1]};

  Point3D Ib   {b[0],    -0.5*b[1],  SQ_3_2*b[1]};
  Point3D IIb  {b[0], -SQ_3_2*b[1],     0.5*b[1]};
  Point3D IIIb {b[0],      -1*b[1],       0*b[1]};
  Point3D IVb  {b[0], -SQ_3_2*b[1],    -0.5*b[1]};
  Point3D Vb   {b[0],    -0.5*b[1], -SQ_3_2*b[1]};
  Point3D VIb  {b[0],       0*b[1],      -1*b[1]};
  Point3D VIIb {b[0],     0.5*b[1], -SQ_3_2*b[1]};
  Point3D VIIIb{b[0],  SQ_3_2*b[1],    -0.5*b[1]};
  Point3D IXb  {b[0],       1*b[1],       0*b[1]};
  Point3D Xb   {b[0],  SQ_3_2*b[1],     0.5*b[1]};
  Point3D XIb  {b[0],     0.5*b[1],  SQ_3_2*b[1]};
  Point3D XIIb {b[0],       0*b[1],       1*b[1]};
  // clang-format on

  // First generation: a single prism
  out[0] = OctType(XIIa, IVb, IVa, VIIIb, VIIIa, XIIb);

  // Second generation: a prism on the outside of each side wall.
  out[1] = OctType(IVb, VIa, VIb, VIIIa, VIIIb, IVa);
  out[2] = OctType(XIIb, IIa, IIb, IVa, IVb, XIIa);
  out[3] = OctType(VIIIb, Xa, Xb, XIIa, XIIb, VIIIa);

  // Third generation: a prism on the outside of each exterior side wall.
  out[4] = OctType(XIIa, Ib, Ia, IIb, IIa, XIIb);
  out[5] = OctType(IIa, IIIb, IIIa, IVb, IVa, IIb);
  out[6] = OctType(IVa, Vb, Va, VIb, VIa, IVb);
  out[7] = OctType(VIa, VIIb, VIIa, VIIIb, VIIIa, VIb);
  out[8] = OctType(VIIIa, IXb, IXa, Xb, Xa, VIIIb);
  out[9] = OctType(Xa, XIb, XIa, XIIb, XIIa, Xb);
}


template<typename ExecPolicy>
void degenerate_segment_test(const char * label, Point2D * polyline, int len, bool expsuccess)
{
  constexpr int gens = 3;
  int octcount = 0;

  SCOPED_TRACE(label);

  OctType * generated = nullptr;
  if (expsuccess)
  {
     EXPECT_TRUE(axom::quest::discretize<ExecPolicy>(polyline, len, gens, generated, octcount));
  }
  else
  {
     EXPECT_FALSE(axom::quest::discretize<ExecPolicy>(polyline, len, gens, generated, octcount));
  }
  EXPECT_EQ(0, octcount);

  axom::deallocate(generated);
}

//------------------------------------------------------------------------------
template<typename ExecPolicy>
void run_degen_segment_tests()
{
  // Test each of the three generations.
  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.

  Point2D * polyline = axom::allocate<Point2D>(2);

  polyline[0] = {0., 0.};
  polyline[1] = {0., 0.};
  degenerate_segment_test<ExecPolicy>("a.x == b.x, a.y == b.y", polyline, 2, true);

  polyline[0] = {1., 0.};
  polyline[1] = {1., 1.};
  degenerate_segment_test<ExecPolicy>("a.x == b.x, a.y != b.y", polyline, 2, true);

  polyline[0] = {1., -0.1};
  polyline[1] = {1.5, 1.};
  degenerate_segment_test<ExecPolicy>("a.y < 0", polyline, 2, false);
  
  polyline[0] = {1., -0.1};
  polyline[1] = {1.5, 1.};
  degenerate_segment_test<ExecPolicy>("a.y < 0", polyline, 2, false);

  polyline[0] = {1., 1.};
  polyline[1] = {1.5, -.1};
  degenerate_segment_test<ExecPolicy>("b.y < 0", polyline, 2, false);

  polyline[0] = {.5, 1.};
  polyline[1] = {0., 1.};
  degenerate_segment_test<ExecPolicy>("a.x > b.x", polyline, 2, false);

  axom::deallocate(polyline);
}

template<typename ExecPolicy>
void one_segment_test(const char * label, Point2D * polyline, int len)
{
  SCOPED_TRACE(label);

  // Test each of the three generations.
  constexpr int generations = 3;

  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int FIRST_GEN_COUNT = 1;
  constexpr int SECOND_GEN_COUNT = 3;
  constexpr int THIRD_GEN_COUNT = 6;

  int generation = 0;

  // The discretized_segment() routine produces a list of 10 hand-calculated
  // octahedra (three generations) that discretize the surface of revolution
  // (SoR) produced by revolving a one-segment polyline around the positive
  // X-axis.
  OctType * handcut = nullptr;
  discretized_segment(polyline[0], polyline[1], handcut);

  OctType * generated = nullptr;
  int octcount = 0;
  axom::quest::discretize<ExecPolicy>(polyline, len, generations, generated, octcount);

  EXPECT_TRUE(
    check_generation(handcut, generated, generation, 0, FIRST_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut,
                               generated,
                               generation,
                               FIRST_GEN_COUNT,
                               SECOND_GEN_COUNT));
  generation += 1;
  EXPECT_TRUE(check_generation(handcut,
                               generated,
                               generation,
                               FIRST_GEN_COUNT + SECOND_GEN_COUNT,
                               THIRD_GEN_COUNT));
  axom::deallocate(generated);
  axom::deallocate(handcut);
}

//------------------------------------------------------------------------------
template<typename ExecPolicy>
void run_single_segment_tests()
{
  Point2D * polyline = axom::allocate<Point2D>(2);

  polyline[0] = Point2D {0.5, 0.};
  polyline[1] = Point2D {1.8, 0.8};
  one_segment_test<ExecPolicy>("Cone (pointing left)", polyline, 2);

  polyline[0] = Point2D {1.0, 1.0};
  polyline[1] = Point2D {1.8, 0.};
  one_segment_test<ExecPolicy>("Cone (pointing right)", polyline, 2);

  polyline[0] = Point2D {1.0, 1.0};
  polyline[1] = Point2D {1.8, 0.8};
  one_segment_test<ExecPolicy>("Truncated cone", polyline, 2);

  polyline[0] = Point2D {1., 1.};
  polyline[1] = Point2D {3., 1.};
  one_segment_test<ExecPolicy>("Cylinder", polyline, 2);

  polyline[0] = Point2D {-.4, 1.2};
  polyline[1] = Point2D {1.2, 1.};
  one_segment_test<ExecPolicy>("a.x < 0, b.x > 0", polyline, 2);

  axom::deallocate(polyline);
}

template<typename ExecPolicy>
void multi_segment_test(const char * label, Point2D * polyline, int len)
{
  SCOPED_TRACE(label);

  // Test each of the three generations.
  constexpr int generations = 3;

  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int FIRST_GEN_COUNT = 1;
  constexpr int SECOND_GEN_COUNT = 3;
  constexpr int THIRD_GEN_COUNT = 6;
  constexpr int TOTAL_COUNT = FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT;

  int generation = 0;

  OctType * generated = nullptr;
  int octcount = 0;
  axom::quest::discretize<ExecPolicy>(polyline, len, generations, generated, octcount);

  int segcount = len-1;
  OctType * handcut = new OctType[segcount * TOTAL_COUNT];

  for (int segidx = 0; segidx < segcount; ++segidx)
  {
     std::stringstream seglabel;
     seglabel << "Segment " << segidx;
     SCOPED_TRACE(seglabel.str());

     // The discretized_segment() routine produces a list of 10 hand-calculated
     // octahedra (three generations) that discretize the surface of revolution
     // (SoR) produced by revolving a one-segment polyline around the positive
     // X-axis.
     OctType * octpointer = &handcut[segidx * TOTAL_COUNT];
     discretized_segment(polyline[segidx], polyline[segidx + 1], octpointer);

     EXPECT_TRUE(check_generation(octpointer,
                                  generated,
                                  generation,
                                  segidx*TOTAL_COUNT + 0,
                                  FIRST_GEN_COUNT));
     generation += 1;
     EXPECT_TRUE(check_generation(octpointer,
                                  generated,
                                  generation,
                                  segidx*TOTAL_COUNT + FIRST_GEN_COUNT,
                                  SECOND_GEN_COUNT));
     generation += 1;
     EXPECT_TRUE(check_generation(octpointer,
                                  generated,
                                  generation,
                                  segidx*TOTAL_COUNT + FIRST_GEN_COUNT + SECOND_GEN_COUNT,
                                  THIRD_GEN_COUNT));
  }

  axom::deallocate(handcut);
  axom::deallocate(generated);
}

//------------------------------------------------------------------------------
template<typename ExecPolicy>
void run_multi_segment_tests()
{
  constexpr int pointcount = 5;
  Point2D * polyline = axom::allocate<Point2D>(pointcount);

  polyline[0] = Point2D {1.0, 0.5};
  polyline[1] = Point2D {1.6, 0.3};
  polyline[2] = Point2D {1.7, 2.0};
  polyline[3] = Point2D {3.1, 1.5};
  polyline[4] = Point2D {4.0, 2.0};
  one_segment_test<ExecPolicy>("multi-segment", polyline, pointcount);

  axom::deallocate(polyline);
}


#endif  /* QUEST_DISCRETIZE_TEST_HPP_ */
