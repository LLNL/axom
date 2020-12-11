// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/inlet_utils.hpp"

namespace axom
{
namespace inlet
{
void setWarningFlag(axom::sidre::Group* root)
{
  if(!root->hasView("warningFlag"))
  {
    root->createViewScalar("warningFlag", 1);
  }
}

std::string appendPrefix(const std::string& prefix, const std::string& name)
{
  return (prefix == "") ? name : prefix + "/" + name;
}

std::string removePrefix(const std::string& prefix, const std::string& name)
{
  if(axom::utilities::string::startsWith(name, prefix + "/"))
  {
    return name.substr(prefix.size());
  }
  SLIC_WARNING(
    fmt::format("[Inlet] Provided name {0} does not "
                "contain prefix {1}",
                name,
                prefix));
  return name;
}

std::string removeBeforeDelimiter(const std::string& path, const char delim)
{
  auto pos = path.find_last_of(delim);
  // Will return an empty string if the delimiter was not found
  return path.substr(pos + 1);
}

std::pair<int, bool> checkedConvertToInt(const std::string& number)
{
  // Use the C versions to avoid the exceptions
  // thrown by std::stoi on conversion failure
  // FIXME: Switch to std::from_chars when C++17 is available
  char* ptr;
  auto result = strtol(number.c_str(), &ptr, 10);
  return {result, *ptr == 0};
}

}  // namespace inlet
}  // namespace axom
