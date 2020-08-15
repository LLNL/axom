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
  // Create intermediate Tables if they don't already exist
  std::string currName = name;
  size_t found = currName.find("/");
  auto currTable = shared_from_this();
  
  while (found != std::string::npos) {
    const std::string& currTableName = getFullName(currTable->m_name, currName.substr(0, found));
    if (!currTable->hasTable(currTableName)) {
      auto table = std::make_shared<Table>(currTableName, "", m_reader, 
                                  m_sidreRootGroup, m_docEnabled);
      currTable->m_tableChildren[currTableName] = table;
      currTable = table;
    } else {
      currTable = currTable->m_tableChildren[currTableName];
    }
    currName = currName.substr(found+1);
    found = currName.find("/");
  }

  std::string currTableName = getFullName(currTable->m_name, currName.substr(0, found));

  if (!currTable->hasTable(currName)) {
    auto table = std::make_shared<Table>(currTableName, description, m_reader, 
                                 m_sidreRootGroup, m_docEnabled);
    currTable->m_tableChildren[currTableName] = table;
    currTable = table;
  } else {
    currTable = currTable->m_tableChildren[currTableName];
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
  std::string fullName = getFullName(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    // This will add any intermediate Tables (if not present) before adding the field
    currTable = addTable(name.substr(0, found));
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
  std::string fullName = getFullName(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    currTable = addTable(name.substr(0, found));
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
  std::string fullName = getFullName(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    currTable = addTable(name.substr(0, found));
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
  std::string fullName = getFullName(m_name, name);
  size_t found =  name.find_last_of("/");

  auto currTable = shared_from_this();
  if (found != std::string::npos) {
    currTable = addTable(name.substr(0, found));
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
  return m_fieldChildren.find(getFullName(m_name, fieldName)) != m_fieldChildren.end();
}

bool Table::hasTable(const std::string& tableName) {
  return m_tableChildren.find(getFullName(m_name, tableName)) != m_tableChildren.end();
}

std::shared_ptr<Table> Table::getTable(std::string name) {
  size_t found = name.find("/");
  auto currTable = shared_from_this();

  while (found != std::string::npos) {
    const std::string& currName = name.substr(0, found);
    if (currTable->hasTable(currName)) {
      currTable = currTable->m_tableChildren[getFullName(currTable->m_name, currName)];
    } else {
      SLIC_WARNING(fmt::format("Table {0} not found", name));
      return nullptr;
    }
    name = name.substr(found+1);
    found = name.find("/");
  }
  
  if (currTable->hasTable(name)) {
    return currTable->m_tableChildren[getFullName(currTable->m_name, name)];
  }
  SLIC_WARNING(fmt::format("Table {0} not found", name));
  return nullptr;
}

std::shared_ptr<Field> Table::getField(std::string name) {
  size_t found =  name.find_last_of("/");
  if (found == std::string::npos) {
    if (hasField(name)) {
      return m_fieldChildren[name];
    }
  } else {
    const std::string& fieldName = name.substr(found+1);
    auto table = getTable(name.substr(0, found));
    if (table && table->hasField(fieldName)) {
      return m_fieldChildren[fieldName];
    }
  }
  SLIC_WARNING(fmt::format("Field {0} not found", name));
  return nullptr;
}

std::string Table::name() {
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
