/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file
 * \brief Defines all-nearest-neighbor queries
 */

#ifndef ANN_QUERY_HPP_
#define ANN_QUERY_HPP_

#include "axom/config.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/UniformGrid.hpp"

namespace axom
{
namespace quest
{

//------------------------------------------------------------------------------
inline double squared_distance(double x1, double y1, double z1,
                               double x2, double y2, double z2)
{
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;

  return dx*dx + dy*dy + dz*dz;
}

//------------------------------------------------------------------------------
void all_nearest_neighbors_bruteforce(double * x, double * y, double * z,
                                      int * region, int n, double limit,
                                      int * neighbor)
{
  // O(n^2) brute force approach.  For each point i, test distance to all other
  // points and report result.

  int * bestsqdist = new int[n];

  for (int i = 0; i < n; ++i) {
    bestsqdist[i] = DBL_MAX;
    neighbor[i] = -1;
  }

  for (int i = 0; i < n; ++i) {
    for (int j = i + 1; j < n; ++j) {
      if (region[i] != region[j]) {
        double sqdist = squared_distance(x[i], y[i], z[i], x[j], y[j], z[j]);
        if (sqdist < bestsqdist[i]) {
          bestsqdist[i] = sqdist;
          neighbor[i] = j;
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
void all_nearest_neighbors_index1(double * x, double * y, double * z,
                                  int * region, int n, double limit,
                                  int * neighbor)
{
  // 1. Build an index, inserting each point individually
  // 2. For each bin A, for each point a in A,
  // 3. For each other bin B, for each point b in B,
  // 4. Compare distances to find the closest distance d = |ab|
  // 5. Exclude a from loop at 2 if we start increasing d, or if D = |AB| > limit.
  // 6. Choose each successive bin B to form shells around A.
}


} // end namespace quest
} // end namespace axom

#endif  // ANN_QUERY_HPP_
