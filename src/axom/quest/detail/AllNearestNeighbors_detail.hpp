/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC.
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

#ifndef ALL_NEAREST_NEIGHBORS_DETAIL_HPP_
#define ALL_NEAREST_NEIGHBORS_DETAIL_HPP_

namespace axom
{
namespace quest
{
namespace detail
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

} // end namespace detail
} // end namespace quest
} // end namespace axom



#endif // ALL_NEAREST_NEIGHBORS_DETAIL_HPP_
