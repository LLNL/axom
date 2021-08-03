// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "quest_discretize.hpp"


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
  for(int i = 0; i < OctType::NUM_OCT_VERTS; ++i)
  {
    Point3D& pt = out[i];
    pt[d] *= -1;
  }
  return out;
}

/* Return a handwritten list of the octahedra discretizing the unit sphere.
 *
 * The routine allocates and returns an array of octahedrons pointed to by out.
 * The caller must free that array.
 */
void discretized_sphere(OctType *& out)
{
  // We're going to return three generations in the out-vector:
  // one in the first generation, eight in the second (covering each
  // of the faces of the first generation), 32 in the third (covering
  // the four exposed faces in each of the second-gen octs).
  constexpr int FIRST_GEN_COUNT = 1;
  constexpr int SECOND_GEN_COUNT = 8;
  constexpr int THIRD_GEN_COUNT = 32;
  out = axom::allocate<OctType>(FIRST_GEN_COUNT + SECOND_GEN_COUNT + THIRD_GEN_COUNT);

  // First generation: one octahedron, with vertices on the unit vectors.
  NAType ihat({1., 0., 0.});
  NAType jhat({0., 1., 0.});
  NAType khat({0., 0., 1.});

  OctType center(Point3D(ihat),
                 Point3D(jhat),
                 Point3D(khat),
                 Point3D(-1 * ihat),
                 Point3D(-1 * jhat),
                 Point3D(-1 * khat));
  out[0] = center;

  // Second generation: eight octs, one for each face of the unit oct.
  // Point ij is halfway between (1, 0, 0) and (0, 1, 0),
  // point jk is halfway between (0, 1, 0) and (0, 0, 1),
  // point ki is halfway between (0, 0, 1) and (1, 0, 0).
  Point3D ij(NAType({M_SQRT1_2, M_SQRT1_2, 0.}));
  Point3D jk(NAType({0., M_SQRT1_2, M_SQRT1_2}));
  Point3D ki(NAType({M_SQRT1_2, 0., M_SQRT1_2}));

  OctType second_gen(Point3D(ihat), Point3D(jhat), Point3D(khat), jk, ki, ij);
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
  double SQRT1_6 = 1. / sqrt(6.);
  // There are three interior points, derived from ij, jk, and ki.
  // Point a is halfway between ij and ki, at (1/sqrt(6))(2, 1, 1).
  Point3D a(SQRT1_6 * NAType({2, 1, 1}));
  // Point b is halfway between ij and jk, at (1/sqrt(6))(1, 2, 1).
  Point3D b(SQRT1_6 * NAType({1, 2, 1}));
  // Point c is halfway between jk and ki, at (1/sqrt(6))(1, 1, 2).
  Point3D c(SQRT1_6 * NAType({1, 1, 2}));

  // There are six edge points, derived from the original corner points and
  // ij, jk, and ki.
  // Point d is halfway between ihat and ij, at
  // (1/sqrt(4 + 2 sqrt(2)))(1+sqrt(2), 1, 0)
  double FACTOR_3G = 1. / sqrt(2. * M_SQRT2 + 4);
  Point3D d(FACTOR_3G * NAType({1. + M_SQRT2, 1., 0}));
  // Point e splits jhat and ij, at (1/sqrt(2 sqrt(2) + 2))(1, 1+sqrt(2), 0)
  Point3D e(FACTOR_3G * NAType({1., 1. + M_SQRT2, 0}));
  // Point f splits jhat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1+sqrt(2), 1)
  Point3D f(FACTOR_3G * NAType({0, 1. + M_SQRT2, 1.}));
  // Point g splits khat and jk, at (1/sqrt(2 sqrt(2) + 2))(0, 1, 1+sqrt(2))
  Point3D g(FACTOR_3G * NAType({0, 1., 1. + M_SQRT2}));
  // Point m splits khat and ki, at (1/sqrt(2 sqrt(2) + 2))(1, 0, 1+sqrt(2))
  Point3D m(FACTOR_3G * NAType({1., 0, 1. + M_SQRT2}));
  // Point n splits ihat and ki, at (1/sqrt(2 sqrt(2) + 2))(1+sqrt(2), 0, 1)
  Point3D n(FACTOR_3G * NAType({1. + M_SQRT2, 0, 1.}));

  int offset = FIRST_GEN_COUNT + SECOND_GEN_COUNT;
  // Here's the first octant of third-generation octahedra (octant 0).
  int octant = 0;
  // First, the interior oct
  out[offset + 0 + octant * 4] = OctType(ij, jk, ki, c, a, b);
  // The one next to ihat
  out[offset + 1 + octant * 4] = OctType(Point3D(ihat), ij, ki, a, n, d);
  // The one next to jhat
  out[offset + 2 + octant * 4] = OctType(Point3D(jhat), jk, ij, b, e, f);
  // The one next to khat
  out[offset + 3 + octant * 4] = OctType(Point3D(khat), ki, jk, c, g, m);

  // Now we get to transform these into all seven remaining octants.
  // Reflect in X
  octant = 1;
  ReflectDimension rd0 = X;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in Y
  octant = 2;
  rd0 = Y;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in Z
  octant = 3;
  rd0 = Z;
  out[offset + 0 + octant * 4] = reflect(rd0, out[offset + 0]);
  out[offset + 1 + octant * 4] = reflect(rd0, out[offset + 1]);
  out[offset + 2 + octant * 4] = reflect(rd0, out[offset + 2]);
  out[offset + 3 + octant * 4] = reflect(rd0, out[offset + 3]);
  // Reflect in X, then Y
  octant = 4;
  rd0 = X;
  ReflectDimension rd1 = Y;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // Reflect in Y, then Z
  octant = 5;
  rd0 = Y;
  rd1 = Z;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // Reflect in Z, then X
  octant = 6;
  rd0 = Z;
  rd1 = X;
  out[offset + 0 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 0]));
  out[offset + 1 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 1]));
  out[offset + 2 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 2]));
  out[offset + 3 + octant * 4] = reflect(rd1, reflect(rd0, out[offset + 3]));
  // And the grand finale: reflect in Z, then Y, then X
  octant = 7;
  rd0 = Z;
  rd1 = Y;
  ReflectDimension rd2 = X;
  out[offset + 0 + octant * 4] =
    reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 0])));
  out[offset + 1 + octant * 4] =
    reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 1])));
  out[offset + 2 + octant * 4] =
    reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 2])));
  out[offset + 3 + octant * 4] =
    reflect(rd2, reflect(rd1, reflect(rd0, out[offset + 3])));
}

//------------------------------------------------------------------------------
TEST(quest_discretize, sphere_test)
{
  // The discretized_sphere() routine produces a list of 41 hand-calculated
  // octahedra (three generations) that discretize the unit sphere.
  OctType * handcut = nullptr;
  discretized_sphere(handcut);

  // The discretize() routine chops up a given sphere into the specified
  // number of generations of octahedra.  Here, we'll discretize the unit
  // sphere in three generations, to match the hand-calculated octs from
  // discretized_sphere().
  SphereType sph;  // Unit sphere at the origin
  constexpr int generations = 3;
  OctType * generated = nullptr;
  int octcount = 0;
  axom::quest::discretize(sph, generations, generated, octcount);

  // Test each of the three generations.
  // We don't know what order they'll be in, but we do know how many octahedra
  // will be in each generation.
  constexpr int FIRST_GEN_COUNT = 1;
  constexpr int SECOND_GEN_COUNT = 8;
  constexpr int THIRD_GEN_COUNT = 32;
  int generation = 0;

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
TEST(quest_discretize, degenerate_sphere_test)
{
  constexpr int generations = 3;
  OctType * generated = nullptr;
  int octcount = 0;

  {
    SCOPED_TRACE("Negative sphere radius");
    SphereType sph(-.2);  // BAD sphere at the origin with negative radius
    EXPECT_FALSE(axom::quest::discretize(sph, generations, generated, octcount));
    EXPECT_EQ(0, octcount);

    axom::deallocate(generated);
  }

  {
    SCOPED_TRACE("Zero sphere radius");
    SphereType sph(0.);  // Degenerate sphere at the origin with zero radius
    EXPECT_TRUE(axom::quest::discretize(sph, generations, generated, octcount));
    EXPECT_EQ(0, octcount);

    axom::deallocate(generated);
  }
}

//------------------------------------------------------------------------------
TEST(quest_discretize, degenerate_segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_degen_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_degen_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && \
  defined(AXOM_USE_UMPIRE) && defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<32> >();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<64> >();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<128> >();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<256> >();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_degen_segment_tests<axom::CUDA_EXEC<512> >();
  }
#endif
}

//------------------------------------------------------------------------------
TEST(quest_discretize, segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_single_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_single_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && \
  defined(AXOM_USE_UMPIRE) && defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<32> >();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<64> >();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<128> >();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<256> >();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_single_segment_tests<axom::CUDA_EXEC<512> >();
  }
#endif
}

//------------------------------------------------------------------------------
TEST(quest_discretize, multi_segment_test)
{
  SLIC_INFO("Discretizing sequentially");
  {
    SCOPED_TRACE("sequential execution");
    run_multi_segment_tests<axom::SEQ_EXEC>();
  }

#if defined(AXOM_USE_OPENMP) && defined(AXOM_USE_RAJA)
  SLIC_INFO("Discretizing with OpenMP");
  {
    SCOPED_TRACE("OpenMP execution");
    run_multi_segment_tests<axom::OMP_EXEC>();
  }
#endif

#if defined(AXOM_USE_CUDA) && defined(AXOM_USE_RAJA) && \
  defined(AXOM_USE_UMPIRE) && defined(__CUDACC__)
  SLIC_INFO("Discretizing with CUDA");
  {
    SCOPED_TRACE("32-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<32> >();
  }

  {
    SCOPED_TRACE("64-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<64> >();
  }

  {
    SCOPED_TRACE("128-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<128> >();
  }

  {
    SCOPED_TRACE("256-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<256> >();
  }

  {
    SCOPED_TRACE("512-wide CUDA execution");
    run_multi_segment_tests<axom::CUDA_EXEC<512> >();
  }
#endif
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
