// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Path.hpp"

#include <iterator>

#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace utilities
{
Path::Path(const std::string& path, const char delim) : m_delim(delim)
{
  // Check if the path has more than one component
  if(path.find(delim) != std::string::npos)
  {
    utilities::string::split(m_components, path, delim);
  }
  else
  {
    m_components.push_back(path);
  }
}

Path Path::join(std::initializer_list<Path> paths, const char delim)
{
  Path result;
  result.m_delim = delim;
  for(const auto& path : paths)
  {
    // Insert everything in order
    result.m_components.insert(result.m_components.end(),
                               std::make_move_iterator(path.m_components.begin()),
                               std::make_move_iterator(path.m_components.end()));
  }
  return result;
}

Path::operator std::string() const
{
  // Start with empty separator to avoid trailing delimiter
  std::string delim_str = "";
  std::string result;
  for(const auto& component : m_components)
  {
    result += delim_str;
    result += component;
    delim_str = m_delim;
  }
  return result;
}

}  // end namespace utilities
}  // end namespace axom
