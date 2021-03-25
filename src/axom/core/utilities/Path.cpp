// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/Path.hpp"

#include <iostream>
#include <iterator>

#include "fmt/fmt.hpp"

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
  else if(!path.empty())
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
    // It seems useful to ignore empty elements here, but maybe this is short-sighted?
    if(!path.m_components.empty())
    {
      result.m_components.insert(
        result.m_components.end(),
        std::make_move_iterator(path.m_components.begin()),
        std::make_move_iterator(path.m_components.end()));
    }
  }
  return result;
}

Path::operator std::string() const
{
  return fmt::format("{0}", fmt::join(m_components, std::string(1, m_delim)));
}

Path Path::parent() const
{
  Path result(*this);
  if(!result.m_components.empty())  // [[likely]] - uncomment when C++20 available
  {
    result.m_components.pop_back();
  }
  return result;
}

std::string Path::baseName() const
{
  if(m_components.empty())
  {
    return {};
  }
  return m_components.back();
}

std::string Path::dirName() const { return static_cast<std::string>(parent()); }

bool operator==(const Path& lhs, const Path& rhs)
{
  return static_cast<std::string>(lhs) == static_cast<std::string>(rhs);
}

}  // end namespace utilities
}  // end namespace axom
