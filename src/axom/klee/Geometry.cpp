// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/Geometry.hpp"
#include "axom/mint/mesh/Mesh.hpp"

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

Geometry::Geometry(const TransformableGeometryProperties &startProperties,
                   axom::sidre::Group* simplexMesh,
                   const std::string& topology,
                   std::shared_ptr<GeometryOperator const> operator_)
  : m_startProperties(startProperties)
  , m_format("memory-blueprint")
  , m_path()
  , m_simplexMesh(simplexMesh)
  , m_topology(topology)
  , m_operator(std::move(operator_))
{
  // axom::mint::Mesh* mintMesh = axom::mint::getMesh(simplexMesh, topo);
  // SLIC_ASSERT(mintMesh->isUnstructured());
  // SLIC_ASSERT(!mintMesh->hasMixedCellTypes());
  // m_simplexMesh.reset(dynamic_cast<SimplexMesh*>(mintMesh));
}

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
