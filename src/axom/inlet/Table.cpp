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


std::shared_ptr<Table> Table::addTable(const std::string& name,
                                       const std::string& description)
{
  std::string fullName = concatenatePaths(m_name, name);
  size_t found = fullName.find("/");
  auto currTable = shared_from_this();
  
  while (found != std::string::npos) {
    const std::string& currName = fullName.substr(0, found);
    if (!currTable->hasTable(currName)) {
      auto table = std::make_shared<Table>(currName, "", m_reader, 
                                  m_sidreRootGroup, m_docEnabled);
      currTable->m_tableChildren[currName] = table;
      currTable = table;
    }
    found = fullName.find("/", found+1);
  }

  if (!currTable->hasTable(fullName)) {
    auto table = std::make_shared<Table>(fullName, description, m_reader, 
                                 m_sidreRootGroup, m_docEnabled);
    currTable->m_tableChildren[fullName] = table;
    currTable = table;
  }

  return currTable;
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
  std::string fullName = concatenatePaths(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    const std::string& tablePath = name.substr(0, found);
    currTable = addTable(tablePath);
  }

  axom::sidre::Group* sidreGroup = currTable->baseFieldAdd(fullName, description);
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

  auto field =  std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::INT8_ID, m_docEnabled);
  currTable->m_fieldChildren[fullName] = field;
  return field;                      
}

std::shared_ptr<Field> Table::addDouble(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = concatenatePaths(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    const std::string& tablePath = name.substr(0, found);
    currTable = addTable(tablePath);
  }

  axom::sidre::Group* sidreGroup = currTable->baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  double value;
  if(m_reader->getDouble(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  auto field = std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::DOUBLE_ID, m_docEnabled);
  currTable->m_fieldChildren[fullName] = field;
  return field;                                   
}

std::shared_ptr<Field> Table::addInt(const std::string& name,
                                     const std::string& description)
{
  std::string fullName = concatenatePaths(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    const std::string& tablePath = name.substr(0, found);
    currTable = addTable(tablePath);
  }

  axom::sidre::Group* sidreGroup = currTable->baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  int value;
  if(m_reader->getInt(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }

  auto field = std::make_shared<axom::inlet::Field>(sidreGroup, m_sidreRootGroup,
                                              axom::sidre::DataTypeId::INT_ID,
                                                                m_docEnabled);
  currTable->m_fieldChildren[fullName] = field;
  return field;
}

std::shared_ptr<Field> Table::addString(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = concatenatePaths(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    const std::string& tablePath = name.substr(0, found);
    currTable = addTable(tablePath);
  }
  axom::sidre::Group* sidreGroup = currTable->baseFieldAdd(fullName, description);
  if (sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }
  
  std::string value;
  if(m_reader->getString(fullName, value))
  {
    sidreGroup->createViewString("value", value);
  }

  auto field = std::make_shared<Field>(sidreGroup, m_sidreRootGroup,
                                 axom::sidre::DataTypeId::CHAR8_STR_ID, m_docEnabled);
  currTable->m_fieldChildren[fullName] = field;
  return field;                                    
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

 std::shared_ptr<Table> Table::registerVerifier(std::function<bool()> lambda) {
  SLIC_WARNING_IF(m_verifier, fmt::format("Verifier for Table {0} already set", 
                                         m_sidreGroup->getPathName()));
  m_verifier = lambda;
  return shared_from_this();
 } 

 bool Table::verify() {
   bool verified = true;
   // Verify this Table
   if (m_verifier && !m_verifier()) {
     verified = false;
     SLIC_WARNING(fmt::format("Table {0} failed verification", m_sidreGroup->getPathName()));
   }
   // Verify the child Fields of this Table
   for (auto field : m_fieldChildren) {
     if (!field.second->verify()) {
       verified = false;
     }
   }
   // Verify the child Tables of this Table
   for (auto table : m_tableChildren) {
     if (!table.second->verify()) {
       verified = false;
     }
   }
   
   return verified;
 }

bool Table::hasField(const std::string& fieldName) {
  return m_fieldChildren.find(concatenatePaths(m_name, fieldName)) != m_fieldChildren.end();
}

bool Table::hasTable(const std::string& tableName) {
  return m_tableChildren.find(concatenatePaths(m_name, tableName)) != m_tableChildren.end();
}

std::string Table::getName() {
  return m_name;
}

std::unordered_map<std::string, std::shared_ptr<Field>> Table::getChildFields() {
  return m_fieldChildren;
}

std::unordered_map<std::string, std::shared_ptr<Table>> Table::getChildTables() {
  return m_tableChildren;
}

} // end namespace inlet
} // end namespace axom
