// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
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
TEST(quest_discretize, sphere_test)
{

   // The discretized_sphere() routine produces a list of 41 hand-calculated
   // octahedra (three generations) that discretize the unit sphere.
   std::vector<OctType> handcut;
   axom::quest::discretized_sphere(handcut);

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
