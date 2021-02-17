// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"

#include "axom/klee/GeometryOperators.hpp"

#include <utility>

namespace axom
{
namespace klee
{
bool operator==(const TransformableGeometryProperties &lhs,
                const TransformableGeometryProperties &rhs)
{
  return lhs.dimensions == rhs.dimensions && lhs.units == rhs.units;
}

Geometry::Geometry(const TransformableGeometryProperties &startProperties,
                   std::string format,
                   std::string path,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format(std::move(format))
  , m_path(std::move(path))
  , m_operator(std::move(operator_))
{ }

TransformableGeometryProperties Geometry::getEndProperties() const
{
  if(m_operator)
  {
    return m_operator->getEndProperties();
  }
  return m_startProperties;
}

}  // namespace klee
}  // namespace axom
