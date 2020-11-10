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

bool checkRequired(const axom::sidre::Group& target, axom::sidre::Group& root)
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

}  // namespace inlet
}  // namespace axom
