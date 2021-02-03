// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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

void setRequired(axom::sidre::Group& target, axom::sidre::Group& root, bool required)
{
  if(target.hasView("required"))
  {
    const std::string msg =
      fmt::format("[Inlet] Required value has already been defined for: {0}",
                  target.getName());

    SLIC_WARNING(msg);
    setWarningFlag(&root);
  }
  else
  {
    if(required)
    {
      target.createViewScalar("required", static_cast<int8>(1));
    }
    else
    {
      target.createViewScalar("required", static_cast<int8>(0));
    }
  }
}

bool checkIfRequired(const axom::sidre::Group& target, axom::sidre::Group& root)
{
  if(!target.hasView("required"))
  {
    return false;
  }
  const axom::sidre::View* valueView = target.getView("required");
  if(valueView == nullptr)
  {
    //TODO: is this possible after it says it has the view?
    return false;
  }
  const int8 intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    const std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      " boolean value named {0}",
      target.getName());
    SLIC_WARNING(msg);
    setWarningFlag(&root);
    return false;
  }

  return static_cast<bool>(intValue);
}

bool verifyRequired(const axom::sidre::Group& target,
                    const bool condition,
                    const std::string& type)
{
  if(target.hasView("required"))
  {
    int8 required = target.getView("required")->getData();
    if(required && !condition)
    {
      const std::string msg = fmt::format(
        "[Inlet] Required {0} not "
        "specified: {1}",
        type,
        target.getPathName());
      SLIC_WARNING(msg);
      return false;
    }
  }
  return true;
}

std::string appendPrefix(const std::string& prefix, const std::string& name)
{
  return (prefix == "") ? name : prefix + "/" + name;
}

std::string removePrefix(const std::string& prefix, const std::string& name)
{
  if(prefix.empty())
  {
    return name;
  }
  else if(axom::utilities::string::startsWith(name, prefix + "/"))
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

std::string removeAllInstances(const std::string& target,
                               const std::string& substr)
{
  std::string result = target;
  auto pos = result.find(substr);
  while(pos != std::string::npos)
  {
    result.erase(pos, substr.length());
    pos = result.find(substr);
  }
  return result;
}

bool checkedConvertToInt(const std::string& number, int& result)
{
  // Use the C versions to avoid the exceptions
  // thrown by std::stoi on conversion failure
  // FIXME: Switch to std::from_chars when C++17 is available
  char* ptr;
  result = strtol(number.c_str(), &ptr, 10);
  return *ptr == 0;
}

void markAsStructContainer(axom::sidre::Group& target)
{
  if(target.hasView(detail::GENERIC_CONTAINER_FLAG))
  {
    // This flag should only ever be one, so we verify that and error otherwise
    const sidre::View* flag = target.getView(detail::GENERIC_CONTAINER_FLAG);
    SLIC_ERROR_IF(
      !flag->isScalar(),
      fmt::format(
        "[Inlet] Generic container flag of group '{0}' was not a scalar",
        target.getName()));
    const int8 value = flag->getScalar();
    SLIC_ERROR_IF(value != 1,
                  fmt::format("[Inlet] Generic container flag of group '{0}' "
                              "had a value other than 1",
                              target.getName()));
  }
  else
  {
    target.createViewScalar(detail::GENERIC_CONTAINER_FLAG, static_cast<int8>(1));
  }
}

}  // namespace inlet
}  // namespace axom
