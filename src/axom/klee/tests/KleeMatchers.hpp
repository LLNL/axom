// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_KLEEMATCHERS_HPP
#define AXOM_KLEEMATCHERS_HPP

#include "gmock/gmock.h"

#include "axom/core/Types.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

namespace axom
{
namespace klee
{
namespace test
{
constexpr double tolerance = 1e-10;

MATCHER_P(AlmostEqMatrix, m, "")
{
  if(arg.getNumRows() != m.getNumRows() ||
     arg.getNumColumns() != m.getNumColumns())
  {
    return false;
  }

  for(IndexType row = 0; row < arg.getNumRows(); ++row)
  {
    for(IndexType column = 0; column < arg.getNumColumns(); ++column)
    {
      if(!::testing::Matches(::testing::DoubleNear(m(row, column), tolerance))(
           arg(row, column)))
      {
        return false;
      }
    }
  }

  return true;
}

MATCHER_P2(AlmostEqArray, array, len, "")
{
  for(int i = 0; i < len; ++i)
  {
    if(!::testing::Matches(::testing::DoubleNear(array[i], tolerance))(arg[i]))
    {
      return false;
    }
  }
  return true;
}

MATCHER_P(AlmostEqVector, vector, "")
{
  for(int i = 0; i < vector.dimension(); ++i)
  {
    if(!::testing::Matches(::testing::DoubleNear(vector[i], tolerance))(arg[i]))
    {
      return false;
    }
  }
  return true;
}

MATCHER_P(AlmostEqPoint, point, "")
{
  for(int i = 0; i < point.dimension(); ++i)
  {
    if(!::testing::Matches(::testing::DoubleNear(point[i], tolerance))(arg[i]))
    {
      return false;
    }
  }
  return true;
}

MATCHER_P(MatchesSlice, slice, "")
{
  return ::testing::Matches(AlmostEqArray(slice.getOrigin(), 3))(arg.getOrigin()) &&
    ::testing::Matches(AlmostEqArray(slice.getNormal(), 3))(arg.getNormal()) &&
    ::testing::Matches(AlmostEqArray(slice.getUp(), 3))(arg.getUp()) &&
    ::testing::Matches(::testing::Eq(slice.getStartProperties()))(
           arg.getStartProperties());
}

}  // namespace test
}  // namespace klee
}  // namespace axom

#endif  //AXOM_KLEEMATCHERS_HPP
