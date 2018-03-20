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

#ifndef ALL_NEAREST_NEIGHBORS_HPP_
#define ALL_NEAREST_NEIGHBORS_HPP_

#include <cfloat>   // for DBL_MAX

#include "axom/config.hpp"

#include "primal/BoundingBox.hpp"
#include "primal/Point.hpp"
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


/*!
 * \brief Find the closest point (in another region) to each given point
 * \param [in] x X-coordinates of input points
 * \param [in] y Y-coordinates of input points
 * \param [in] z Z-coordinates of input points
 * \param [in] region Region of each point
 * \param [in] n Number of points
 * \param [in] limit Max distance for all-nearest-neighbors query
 * \param [out] neighbor Index of nearest neighbor not in the same class
 * \param [out] sqdistance Squared distance to nearest neighbor
 * \pre x, y, z, and region have n entries
 * \pre neighbor is allocated with room for n entries
 *
 * This method compares each point p to each other point in the UniformGrid
 * bins that are at least partly within a Manhattan distance of limit
 * from p.  This cuts out the far-away points.  We expect this will result
 * in a substantial time savings over the brute-force algorithm, but the
 * run time is dependent on the point distribution.
 */
void all_nearest_neighbors(const double* x, const double* y, const double* z,
                           const int* region, int n, double limit,
                           int* neighbor, double* sqdistance)
{
  // Indexed approach.  For each point i, test distance to all other
  // points in this and neighboring UniformGrid bins (out to distance limit)
  // and report result.

  typedef axom::primal::UniformGrid<int, 3> GridType;
  typedef GridType::BoxType BoxType;
  typedef GridType::PointType PointType;
  typedef BoxType::VectorType VectorType;

  double sqlimit = limit * limit;

  PointType pmin, pmax;
  for (int i = 0 ; i < n ; ++i)
  {
    sqdistance[i] = DBL_MAX;
    neighbor[i] = -1;
    pmin[0] = std::min(pmin[0], x[i]);
    pmax[0] = std::max(pmax[0], x[i]);
    pmin[1] = std::min(pmin[1], y[i]);
    pmax[1] = std::max(pmax[1], y[i]);
    pmin[2] = std::min(pmin[2], z[i]);
    pmax[2] = std::max(pmax[2], z[i]);
  }
  BoxType allpointsbox(pmin, pmax);

  int res[3];
  VectorType boxrange = allpointsbox.range();
  for (int i = 0 ; i < 3 ; ++i)
  {
    res[i] = std::max(1, (int)(boxrange[i] / limit + 0.5));
  }

  // 1. Build an index, inserting each point individually
  GridType ugrid(allpointsbox, res);
  for (int i = 0 ; i < n ; ++i)
  {
    ugrid.insert(BoxType(PointType::make_point(x[i], y[i], z[i])), i);
  }

  // 2. For point a,
  for (int i = 0 ; i < n ; ++i)
  {
    // 3. For each other bin B less than limit distance away from a, for each
    // point b in B,
    PointType qmin =
      PointType::make_point(x[i] - limit, y[i] - limit, z[i] - limit);
    PointType qmax =
      PointType::make_point(x[i] + limit, y[i] + limit, z[i] + limit);
    BoxType qbox(qmin, qmax);
    const std::vector<int> qbins = ugrid.getBinsForBbox(qbox);
    const size_t querybincount = qbins.size();
    for (size_t binidx = 0 ; binidx < querybincount ; ++binidx)
    {
      const std::vector<int> bs = ugrid.getBinContents(qbins[binidx]);
      const size_t binsize = bs.size();
      for (size_t bj = 0 ; bj < binsize ; ++bj)
      {
        // 4. Compare distances to find the closest distance d = |ab|
        int j = bs[bj];
        if (region[i] != region[j])
        {
          double sqdist = squared_distance(x[i], y[i], z[i], x[j], y[j], z[j]);
          if (sqdist < sqdistance[i] && sqdist < sqlimit)
          {
            sqdistance[i] = sqdist;
            neighbor[i] = j;
          }
        }
      }
    }
  }
}


} // end namespace quest
} // end namespace axom

#endif  // ALL_NEAREST_NEIGHBORS_HPP_
