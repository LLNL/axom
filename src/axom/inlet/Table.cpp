// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Table.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"
#include "axom/inlet/inlet_utils.hpp"

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

std::shared_ptr<Table> Table::addTable(const std::string& name,
                                       const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  return std::make_shared<Table>(fullName, description, m_reader, 
                                 m_sidreRootGroup, m_docEnabled);
}

axom::sidre::Group* Table::baseFieldAdd(const std::string& name,
                                        const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "Inlet's Reader class not set");
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  if (m_sidreRootGroup->hasGroup(name))
  {
    SLIC_WARNING("Inlet: Cannot add value that already exists: " + name);
    setWarningFlag(m_sidreRootGroup);
    return nullptr;
  }

  axom::sidre::Group* sidreGroup = m_sidreRootGroup->createGroup(name);
  SLIC_ASSERT_MSG(sidreGroup != nullptr, "Sidre failed to create group");
  if (description != "")
  {
    sidreGroup->createViewString("description", description);
  }

  return sidreGroup;
}

std::shared_ptr<Field> Table::addBool(const std::string& name,
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
    sidreGroup->createViewScalar("value", value? int8(1) : int8(0) );
  }

  return std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::INT8_ID, m_docEnabled);
}

std::shared_ptr<Field> Table::addDouble(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  double value;
  if(m_reader->getDouble(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  return std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::DOUBLE_ID, m_docEnabled);
}

std::shared_ptr<Field> Table::addInt(const std::string& name,
                                     const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  int value;
  if(m_reader->getInt(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  return std::make_shared<axom::inlet::Field>(sidreGroup, m_sidreRootGroup,
                                              axom::sidre::DataTypeId::INT_ID,
                                                                m_docEnabled);
}

std::shared_ptr<Field> Table::addString(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = getFullName(m_name, name);
  axom::sidre::Group* sidreGroup = baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  std::string value;
  if(m_reader->getString(fullName, value))
  {
    sidreGroup->createViewString("value", value);
  }

  return std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::CHAR8_STR_ID, m_docEnabled);
}

std::shared_ptr<Table> Table::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Table specific Sidre Datastore Group not set");

  if (m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format("Inlet Table has already defined required value: {0}",
                                  m_sidreGroup->getName());
    
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return shared_from_this();
  }

  if (isRequired)
  {
    m_sidreGroup->createViewScalar("required", (int8)1);
  }
  else
  {
    m_sidreGroup->createViewScalar("required", (int8)0);
  }

  return shared_from_this();
}

bool Table::required()
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Table specific Sidre Datastore Group not set");

  if (!m_sidreGroup->hasView("required"))
  {
    return false;
  }
  axom::sidre::View* valueView = m_sidreGroup->getView("required");
  if (valueView == nullptr)
  {
    //TODO: is this possible after it says it has the view?
    return false;
  }
  int8 intValue = valueView->getScalar();
  if (intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format("Invalid integer value stored in boolean"
                                  " value named {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return false;
  }

  return (bool)intValue;
}

} // end namespace inlet
} // end namespace axom
