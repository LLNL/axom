// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Group.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

std::string getFullName(const std::string& prefix, const std::string& name)
{
  if (prefix == "")
  {
    return name;
  }
  else
  {
    return prefix + "/" + name;
  }
}

std::shared_ptr<Group> Group::addGroup(const std::string& name,
                                       const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  return std::make_shared<Group>(fullName, description, m_reader, m_sidreRootGroup);
}

axom::sidre::Group* Group::baseFieldAdd(const std::string& name,
                                        const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  if (m_sidreRootGroup->hasGroup(name))
  {
    SLIC_WARNING("Inlet: Cannot add value that already exists: " + name);
    return nullptr;
  }

  axom::sidre::Group* sidreGroup = m_sidreRootGroup->createGroup(name);
  SLIC_ASSERT_MSG(sidreGroup != nullptr, "Sidre failed to create group");
  if (description == "")
  {
    sidreGroup->createViewString("description", description);
  }

  return sidreGroup;
}

std::shared_ptr<Field> Group::addBool(const std::string& name,
                                      const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    //TODO: better idea?
    return std::shared_ptr<Field>(nullptr);
  }

  bool value;
  if(m_reader->getBool(fullName, value))
  {
    if (value)
    {
      sidreGroup->createViewScalar("value", (int8)1);
    }
    else
    {
      sidreGroup->createViewScalar("value", (int8)0);
    }
  }

  return std::make_shared<Field>(sidreGroup);
}

std::shared_ptr<Field> Group::addDouble(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  double value;
  if(m_reader->getDouble(name, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  return std::make_shared<Field>(sidreGroup);
}

std::shared_ptr<Field> Group::addInt(const std::string& name,
                                     const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  int value;
  if(m_reader->getInt(name, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  return std::make_shared<axom::inlet::Field>(sidreGroup);
}

std::shared_ptr<Field> Group::addString(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  std::string value;
  if(m_reader->getString(name, value))
  {
    sidreGroup->createViewString("value", value);
  }

  return std::make_shared<Field>(sidreGroup);
}

} // end namespace inlet
} // end namespace axom
