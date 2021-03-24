// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/ShapeSet.hpp"

#include <utility>
#include <stdexcept>

#include "axom/core/utilities/FileUtilities.hpp"

namespace axom
{
namespace klee
{
void ShapeSet::setShapes(std::vector<Shape> shapes)
{
  m_shapes = std::move(shapes);
}

void ShapeSet::setPath(const std::string &path) { m_path = path; }

std::string ShapeSet::resolvePath(const std::string &filePath) const
{
  if(m_path.empty())
  {
    throw std::logic_error("The ShapeSet's path has not been set");
  }
  if(filePath[0] == '/')
  {
    return filePath;
  }
  std::string dir;
  utilities::filesystem::getDirName(dir, m_path);
  return utilities::filesystem::joinPath(dir, filePath);
}
}  // namespace klee
}  // namespace axom
