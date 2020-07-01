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

std::shared_ptr<Group> Inlet::addGroup(const std::string& name,
                                       const std::string& description)
{
  return std::make_shared<Group>(name, description, m_reader, m_sidreRootGroup);
}

std::shared_ptr<Field> Inlet::addBool(const std::string& name,
                                      const std::string& description)
{
  return m_group->addBool(name, description);
}

std::shared_ptr<Field> Inlet::addDouble(const std::string& name,
                                        const std::string& description)
{
  return m_group->addDouble(name, description);
}

std::shared_ptr<Field> Inlet::addInt(const std::string& name,
                                     const std::string& description)
{
  return m_group->addInt(name, description);
}

std::shared_ptr<Field> Inlet::addString(const std::string& name,
                                        const std::string& description)
{
  return m_group->addString(name, description);
}


//-------------------------------------------------
//   Get values out of the datastore
//-------------------------------------------------

axom::sidre::View* Inlet::baseGet(const std::string& name)
{
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  // All data hangs under the group's name
  if (!m_sidreRootGroup->hasGroup(name))
  {
    return nullptr;
  }
  axom::sidre::Group* group = m_sidreRootGroup->getGroup(name);
  if (group == nullptr)
  {
    return nullptr;
  }

  if (!group->hasView("value"))
  {
    return nullptr;
  }
  return group->getView("value");
}

bool Inlet::get(const std::string& name, bool& value)
{
  axom::sidre::View* valueView = baseGet(name);
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
                                  name);
    SLIC_WARNING(msg);
    return false;
  }

  value = (bool)intValue;
  return true;
}

bool Inlet::get(const std::string& name, double& value)
{
  axom::sidre::View* valueView = baseGet(name);
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
  axom::sidre::View* valueView = baseGet(name);
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
  axom::sidre::View* valueView = baseGet(name);
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
