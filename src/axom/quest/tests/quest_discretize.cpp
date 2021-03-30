// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/slic/interface/slic.hpp"

#include "axom/quest/Discretize.hpp"  // quest::Discretize

// Google Test includes
#include "gtest/gtest.h"

#include <vector>

using SphereType = axom::primal::Sphere<double, 3>;
using OctType = axom::primal::Octahedron<double, 3>;
using PointType = axom::primal::Point<double, 3>;
using TwoDPointType = axom::primal::Point<double, 2>;
using NAType = axom::primal::NumericArray<double, 3>;

//------------------------------------------------------------------------------
bool check_generation(std::vector<OctType> & standard,
                      std::vector<OctType> & test,
                      int generation,
                      int offset,
                      int count)
{
   int *matched = new int[count];
   // Clear the array of "matched" flags
   for (int i = 0; i < count; ++i) { matched[i] = 0; }

   for (int i = 0; i < count; ++i)
   {
      int has_matched = -1;
      for (int j = 0; j < count; ++j)
      {
         if (!matched[j] && standard[offset+i].equals(test[offset+j]))
         {
            matched[j] = 1;
            has_matched = offset+j;
         }
      }
      if (has_matched < 0)
      {
         std::cout << "Gen " << generation << " standard oct " <<
            (offset+i) << " didn't match: " << standard[offset+i] <<
            std::endl;
      }
   }

   int matchcount = 0;
   for (int i = 0; i < count; ++i) { matchcount += matched[i]; }

   bool matches = (matchcount == count);

   if (!matches)
   {
      // a little more reporting here
      std::cout << "Generation " << generation << " had " <<
         (count - matchcount) <<
         " test octahedra not matched to standard octahedra:" << std::endl;
      for (int i = 0; i < count; ++i)
      {
         if (!matched[i])
         {
            std::cout << "Test oct " << (offset + i) << " not matched:" <<
               test[offset+i] << std::endl;
         }
      }
   }

   delete[] matched;
   return matches;
}

//------------------------------------------------------------------------------
enum ReflectDimension
{
   X = 0,
   Y,
   Z
};

/* Reflect an octahedron in a given dimension.
 */
OctType reflect(ReflectDimension d, OctType o)
{
   OctType out(o);
   for (int i = 0; i < OctType::NUM_OCT_VERTS; ++i)
   {
      PointType & pt = out[i];
      pt[d] *= -1;
   }
   return out;
}

/* Return a handwritten list of the octahedra discretizing the unit sphere.
 */
void discretized_sphere(std::vector<OctType> & out)
{
   // We're going to return three generations in the out-vector:
   // one in the first generation, eight in the second (covering each
   // of the faces of the first generation), 32 in the third (covering
   // the four exposed faces in each of the second-gen octs).
   constexpr int FIRST_GEN_COUNT = 1;
   constexpr int SECOND_GEN_COUNT = 8;
   constexpr int THIRD_GEN_COUNT = 32;
   out.resize(FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT);
   
   // First generation: one octahedron, with vertices on the unit vectors.
   NAType ihat({1., 0., 0.});
   NAType jhat({0., 1., 0.});
   NAType khat({0., 0., 1.});

   OctType center(PointType(ihat), PointType(jhat), PointType(khat),
                  PointType(-1 * ihat), PointType(-1 * jhat), 
                  PointType(-1 * khat));
   out[0] = center;
   
   // Second generation: eight octs, one for each face of the unit oct.
   // Point ij is halfway between (1, 0, 0) and (0, 1, 0),
   // point jk is halfway between (0, 1, 0) and (0, 0, 1),
   // point ki is halfway between (0, 0, 1) and (1, 0, 0).
   PointType ij(NAType({M_SQRT1_2, M_SQRT1_2, 0.}));
   PointType jk(NAType({0., M_SQRT1_2, M_SQRT1_2}));
   PointType ki(NAType({M_SQRT1_2, 0., M_SQRT1_2}));

   OctType second_gen(PointType(ihat), PointType(jhat), PointType(khat),
                      jk, ki, ij);
   out[1] = second_gen;
   out[2] = reflect(X, second_gen);
   out[3] = reflect(Y, second_gen);
   out[4] = reflect(Z, second_gen);
   out[5] = reflect(Y, reflect(X, second_gen));
   out[6] = reflect(Z, reflect(Y, second_gen));
   out[7] = reflect(X, reflect(Z, second_gen));
   out[8] = reflect(X, reflect(Y, reflect(Z, second_gen)));
   
   // Third generation: 32 new octs, one for each exposed face of the previous
   // generation.
   double SQRT1_6 = 1./sqrt(6.);
   // There are three interior points, derived from ij, jk, and ki.
   // Point a is halfway between ij and ki, at (1/sqrt(6))(2, 1, 1).
   PointType a(SQRT1_6*NAType({2, 1, 1}));
   // Point b is halfway between ij and jk, at (1/sqrt(6))(1, 2, 1).
   PointType b(SQRT1_6*NAType({1, 2, 1}));
   // Point c is halfway between jk and ki, at (1/sqrt(6))(1, 1, 2).
   PointType c(SQRT1_6*NAType({1, 1, 2}));
   
   // There are six edge points, derived from the original corner points and
   // ij, jk, and ki.
   // Point d is halfway between ihat and ij, at
   // (1/sqrt(4 + 2 sqrt(2)))(1+sqrt(2), 1, 0)
   double FACTOR_3G = 1./sqrt(2.*M_SQRT2 + 4);
   PointType d(FACTOR_3G*NAType({1. + M_SQRT2, 1., 0}));
   // Point e splits jhat and ij, at (1/sqrt(2 sqrt(2) + 2))(1, 1+sqrt(2), 0)
   PointType e(FACTOR_3G*NAType({1., 1. + M_SQRT2, 0}));
   // Point f splits jhat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1+sqrt(2), 1)
   PointType f(FACTOR_3G*NAType({0, 1. + M_SQRT2, 1.}));
   // Point g splits khat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1, 1+sqrt(2))
   PointType g(FACTOR_3G*NAType({0, 1., 1. + M_SQRT2}));
   // Point m splits khat and ki, at (1/sqrt(2 sqrt(2) + 2))(1, 0, 1+sqrt(2))
   PointType m(FACTOR_3G*NAType({1., 0, 1. + M_SQRT2}));
   // Point n splits ihat and ki, at (1/sqrt(2 sqrt(2) + 2))(1+sqrt(2), 0, 1)
   PointType n(FACTOR_3G*NAType({1. + M_SQRT2, 0, 1.}));

   int offset = FIRST_GEN_COUNT + SECOND_GEN_COUNT;
   // Here's the first octant of third-generation octahedra (octant 0).
   int octant = 0;
   // First, the interior oct
   out[offset+0+octant*4] = OctType(ij, jk, ki, c, a, b);
   // The one next to ihat
   out[offset+1+octant*4] = OctType(PointType(ihat), ij, ki, a, n, d);
   // The one next to jhat
   out[offset+2+octant*4] = OctType(PointType(jhat), jk, ij, b, e, f);
   // The one next to khat
   out[offset+3+octant*4] = OctType(PointType(khat), ki, jk, c, g, m);

   // Now we get to transform these into all seven remaining octants.
   // Reflect in X
   octant = 1;
   ReflectDimension rd0 = X;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in Y
   octant = 2;
   rd0 = Y;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in Z
   octant = 3;
   rd0 = Z;
   out[offset+0+octant*4] = reflect(rd0, out[offset+0]);
   out[offset+1+octant*4] = reflect(rd0, out[offset+1]);
   out[offset+2+octant*4] = reflect(rd0, out[offset+2]);
   out[offset+3+octant*4] = reflect(rd0, out[offset+3]);
   // Reflect in X, then Y
   octant = 4;
   rd0 = X;
   ReflectDimension rd1 = Y;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // Reflect in Y, then Z
   octant = 5;
   rd0 = Y;
   rd1 = Z;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // Reflect in Z, then X
   octant = 6;
   rd0 = Z;
   rd1 = X;
   out[offset+0+octant*4] = reflect(rd1, reflect(rd0, out[offset+0]));
   out[offset+1+octant*4] = reflect(rd1, reflect(rd0, out[offset+1]));
   out[offset+2+octant*4] = reflect(rd1, reflect(rd0, out[offset+2]));
   out[offset+3+octant*4] = reflect(rd1, reflect(rd0, out[offset+3]));
   // And the grand finale: reflect in Z, then Y, then X
   octant = 7;
   rd0 = Z;
   rd1 = Y;
   ReflectDimension rd2 = X;
   out[offset+0+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+0])));
   out[offset+1+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+1])));
   out[offset+2+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+2])));
   out[offset+3+octant*4] = reflect(rd2, reflect(rd1, reflect(rd0, out[offset+3])));
}

/* Return a handwritten list of the octahedra discretizing a one-segment polyline.
 */
void discretized_segment(TwoDPointType a, TwoDPointType b, std::vector<OctType> & out)
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
   out.resize(FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT);

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
   
   const double SQ_3_2 = sqrt(3.)/2.;

   PointType Ia   {a[0],    -0.5*a[1],  SQ_3_2*a[1]};
   PointType IIa  {a[0], -SQ_3_2*a[1],     0.5*a[1]};
   PointType IIIa {a[0],      -1*a[1],       0*a[1]};
   PointType IVa  {a[0], -SQ_3_2*a[1],    -0.5*a[1]};
   PointType Va   {a[0],    -0.5*a[1], -SQ_3_2*a[1]};
   PointType VIa  {a[0],       0*a[1],      -1*a[1]};
   PointType VIIa {a[0],     0.5*a[1], -SQ_3_2*a[1]};
   PointType VIIIa{a[0],  SQ_3_2*a[1],    -0.5*a[1]};
   PointType IXa  {a[0],       1*a[1],       0*a[1]};
   PointType Xa   {a[0],  SQ_3_2*a[1],     0.5*a[1]};
   PointType XIa  {a[0],     0.5*a[1],  SQ_3_2*a[1]};
   PointType XIIa {a[0],       0*a[1],       1*a[1]};

   PointType Ib   {b[0],    -0.5*b[1],  SQ_3_2*b[1]};
   PointType IIb  {b[0], -SQ_3_2*b[1],     0.5*b[1]};
   PointType IIIb {b[0],      -1*b[1],       0*b[1]};
   PointType IVb  {b[0], -SQ_3_2*b[1],    -0.5*b[1]};
   PointType Vb   {b[0],    -0.5*b[1], -SQ_3_2*b[1]};
   PointType VIb  {b[0],       0*b[1],      -1*b[1]};
   PointType VIIb {b[0],     0.5*b[1], -SQ_3_2*b[1]};
   PointType VIIIb{b[0],  SQ_3_2*b[1],    -0.5*b[1]};
   PointType IXb  {b[0],       1*b[1],       0*b[1]};
   PointType Xb   {b[0],  SQ_3_2*b[1],     0.5*b[1]};
   PointType XIb  {b[0],     0.5*b[1],  SQ_3_2*b[1]};
   PointType XIIb {b[0],       0*b[1],       1*b[1]};

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

//------------------------------------------------------------------------------
TEST(quest_discretize, sphere_test)
{

   // The discretized_sphere() routine produces a list of 41 hand-calculated
   // octahedra (three generations) that discretize the unit sphere.
   std::vector<OctType> handcut;
   discretized_sphere(handcut);

   // The discretize() routine chops up a given sphere into the specified
   // number of generations of octahedra.  Here, we'll discretize the unit
   // sphere in three generations, to match the hand-calculated octs from
   // discretized_sphere().
   SphereType sph; // Unit sphere at the origin
   constexpr int generations = 3;
   std::vector<OctType> generated;
   axom::quest::discretize(sph, generations, generated);

   // Test each of the three generations.
   // We don't know what order they'll be in, but we do know how many octahedra
   // will be in each generation.
   constexpr int FIRST_GEN_COUNT = 1;
   constexpr int SECOND_GEN_COUNT = 8;
   constexpr int THIRD_GEN_COUNT = 32;
   int generation = 0;

   EXPECT_TRUE(check_generation(handcut,
                                generated,
                                generation,
                                0,
                                FIRST_GEN_COUNT));
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
}

//------------------------------------------------------------------------------
TEST(quest_discretize, segment_test)
{

   // The discretize() routine chops up a given SoR into the specified
   // number of generations of octahedra.  Here, we'll discretize the polyline
   // in the previous comment in three generations, to match the hand-calculated
   // octs from discretized_segment().
   TwoDPointType a{1.0, 1.0};
   TwoDPointType b{1.8, 0.8};
   std::vector<TwoDPointType> polyline;
   polyline.push_back(a);
   polyline.push_back(b);

   // The discretized_segment() routine produces a list of 10 hand-calculated
   // octahedra (three generations) that discretize the surface of revolution
   // (SoR) produced by revolving a one-segment polyline around the positive
   // X-axis.
   std::vector<OctType> handcut;
   discretized_segment(a, b, handcut);

   constexpr int generations = 3;
   std::vector<OctType> generated;
   axom::quest::discretize(polyline, generations, generated);

   // Test each of the three generations.
   // We don't know what order they'll be in, but we do know how many octahedra
   // will be in each generation.
   constexpr int FIRST_GEN_COUNT = 1;
   constexpr int SECOND_GEN_COUNT = 3;
   constexpr int THIRD_GEN_COUNT = 6;
   
   int generation = 0;
   EXPECT_TRUE(check_generation(handcut,
                                generated,
                                generation,
                                0,
                                FIRST_GEN_COUNT));
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
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
