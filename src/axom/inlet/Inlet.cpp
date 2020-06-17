// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

axom::sidre::Group* Inlet::addGroup(const std::string& name,
                                    const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Sidre Datastore class not set");

  if (m_sidreGroup->hasGroup(name))
  {
    SLIC_WARNING("Inlet: Cannot create Group that already exists: " + name);
    return nullptr;
  }
  axom::sidre::Group* group = m_sidreGroup->createGroup(name);
  group->createViewString("description", description);

  return group;
}

axom::sidre::Group* Inlet::add(const std::string& name,
                               const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  if (m_sidreGroup->hasGroup(name))
  {
    SLIC_WARNING("Inlet: Cannot add value that already exists: " + name);
    return nullptr;
  }

  axom::sidre::Group* group = m_sidreGroup->createGroup(name);
  SLIC_ASSERT_MSG(group != nullptr, "Sidre failed to create group");
  if (description == "")
  {
    group->createViewString("description", description);
  }

  return group;
}

axom::sidre::Group* Inlet::addBool(const std::string& name,
                                   const std::string& description)
{
  axom::sidre::Group* group = add(name, description);
  if (group == nullptr)
  {
    return nullptr;
  }

  bool value;
  if(m_reader->getBool(name, value))
  {
    if (value)
    {
      group->createViewScalar("value", (int8)1);
    }
    else
    {
      group->createViewScalar("value", (int8)0);
    }
  }

  return group;
}

axom::sidre::Group* Inlet::addDouble(const std::string& name,
                                     const std::string& description)
{
  axom::sidre::Group* group = add(name, description);
  if (group == nullptr)
  {
    return nullptr;
  }
  
  double value;
  if(m_reader->getDouble(name, value))
  {
    group->createViewScalar("value", value);
  }

  return group;
}

axom::sidre::Group* Inlet::addInt(const std::string& name,
                                  const std::string& description)
{
  axom::sidre::Group* group = add(name, description);
  if (group == nullptr)
  {
    return nullptr;
  }
  
  int value;
  if(m_reader->getInt(name, value))
  {
    group->createViewScalar("value", value);
  }

  return group;
}

axom::sidre::Group* Inlet::addString(const std::string& name,
                                     const std::string& description)
{
  axom::sidre::Group* group = add(name, description);
  if (group == nullptr)
  {
    return nullptr;
  }
  
  std::string value;
  if(m_reader->getString(name, value))
  {
    group->createViewString("value", value);
  }

  return group;
}


//-------------------------------------------------
//   Get values out of the datastore
//-------------------------------------------------

axom::sidre::View* Inlet::get(const std::string& name)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  // All data hangs under the group's name
  if (!m_sidreGroup->hasGroup(name))
  {
    return nullptr;
  }
  axom::sidre::Group* group = m_sidreGroup->getGroup(name);
  if (group == nullptr)
  {
    return nullptr;
  }

  if (!group->hasView("value"))
  {
    return nullptr;
  }
  axom::sidre::View* valueView = group->getView("value");
  if (valueView == nullptr)
  {
    return nullptr;
  }
  return valueView;
}

bool Inlet::get(const std::string& name, bool& value)
{
  axom::sidre::View* valueView = get(name);
  if (valueView == nullptr)
  {
    return false;
  }

  // There is no boolean type in conduit/sidre so we use int8
  if (valueView->getTypeID() != axom::sidre::INT8_ID)
  {
    std::string msg = fmt::format("Boolean named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    SLIC_WARNING(msg);
    return false;
  }

  int8 intValue = valueView->getScalar();
  if (intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format("Invalid integer value stored in boolean"
                                  " value named {0}",
                                  name, valueView->getTypeID());
    SLIC_WARNING(msg);
    return false;
  }

  value = (bool)intValue;
  return true;
}

bool Inlet::get(const std::string& name, double& value)
{
  axom::sidre::View* valueView = get(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::DOUBLE_ID)
  {
    std::string msg = fmt::format("Double named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    SLIC_WARNING(msg);
    return false;
  }

  value = valueView->getScalar();
  return true;
}

bool Inlet::get(const std::string& name, int& value)
{
  axom::sidre::View* valueView = get(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::INT_ID)
  {
    std::string msg = fmt::format("Integer named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    SLIC_WARNING(msg);
    return false;
  }

  value = valueView->getScalar();
  return true;
}

bool Inlet::get(const std::string& name, std::string& value)
{
  axom::sidre::View* valueView = get(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::CHAR8_STR_ID)
  {
    std::string msg = fmt::format("String named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    SLIC_WARNING(msg);
    return false;
  }

  const char* valueStr = valueView->getString();
  if (valueStr == nullptr)
  {
    value = std::string("");
  }
  value = std::string(valueStr);
  return true;
}

} // end namespace inlet
} // end namespace axom
