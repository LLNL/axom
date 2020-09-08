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

  while(found != std::string::npos)
  {
    const std::string& currTableName =
      appendPrefix(currTable->m_name, currName.substr(0, found));
    if(!currTable->hasChildTable(currTableName))
    {
      auto table = std::make_shared<Table>(currTableName,
                                           "",
                                           m_reader,
                                           m_sidreRootGroup,
                                           m_docEnabled);
      currTable->m_tableChildren[currTableName] = table;
      currTable = table;
    }
    else
    {
      currTable = currTable->m_tableChildren[currTableName];
    }
    currName = currName.substr(found + 1);
    found = currName.find("/");
  }

  std::string currTableName =
    appendPrefix(currTable->m_name, currName.substr(0, found));

  if(!currTable->hasChildTable(currName))
  {
    auto table = std::make_shared<Table>(currTableName,
                                         description,
                                         m_reader,
                                         m_sidreRootGroup,
                                         m_docEnabled);
    currTable->m_tableChildren[currTableName] = table;
    currTable = table;
  }
  else
  {
    currTable = currTable->m_tableChildren[currTableName];
  }

  return currTable;
}

axom::sidre::Group* Table::createSidreGroup(const std::string& name,
                                            const std::string& description)
{
  SLIC_ASSERT_MSG(m_reader != nullptr, "[Inlet] Reader class not set");
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr,
                  "[Inlet] Sidre Datastore Group not set");

  if(m_sidreRootGroup->hasGroup(name))
  {
    SLIC_WARNING("[Inlet] Cannot add value that already exists: " + name);
    setWarningFlag(m_sidreRootGroup);
    return nullptr;
  }

  axom::sidre::Group* sidreGroup = m_sidreRootGroup->createGroup(name);
  sidreGroup->createViewString("InletType", "Field");
  SLIC_ASSERT_MSG(sidreGroup != nullptr, "[Inlet] Sidre failed to create group");
  if(description != "")
  {
    sidreGroup->createViewString("description", description);
  }

  return sidreGroup;
}

std::shared_ptr<Field> Table::addField(axom::sidre::Group* sidreGroup,
                                       axom::sidre::DataTypeId type,
                                       const std::string& fullName,
                                       const std::string& name)
{
  size_t found = name.find_last_of("/");
  auto currTable = shared_from_this();
  if(found != std::string::npos)
  {
    // This will add any intermediate Tables (if not present) before adding the field
    currTable = addTable(name.substr(0, found));
  }
  auto field = std::make_shared<axom::inlet::Field>(sidreGroup,
                                                    m_sidreRootGroup,
                                                    type,
                                                    m_docEnabled);
  currTable->m_fieldChildren[fullName] = field;
  return field;
}

std::shared_ptr<Field> Table::addBool(const std::string& name,
                                      const std::string& description)
{
  std::string fullName = appendPrefix(m_name, name);
  axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
  if(sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }

  bool value;
  if(m_reader->getBool(fullName, value))
  {
    sidreGroup->createViewScalar("value", value ? int8(1) : int8(0));
  }
  return addField(sidreGroup, axom::sidre::DataTypeId::INT8_ID, fullName, name);
}

std::shared_ptr<Field> Table::addDouble(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = appendPrefix(m_name, name);
  axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
  if(sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }

  double value;
  if(m_reader->getDouble(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }
  return addField(sidreGroup, axom::sidre::DataTypeId::DOUBLE_ID, fullName, name);
}

std::shared_ptr<Field> Table::addInt(const std::string& name,
                                     const std::string& description)
{
  std::string fullName = appendPrefix(m_name, name);
  axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
  if(sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }

  int value;
  if(m_reader->getInt(fullName, value))
  {
    sidreGroup->createViewScalar("value", value);
  }
  return addField(sidreGroup, axom::sidre::DataTypeId::INT_ID, fullName, name);
}

std::shared_ptr<Field> Table::addString(const std::string& name,
                                        const std::string& description)
{
  std::string fullName = appendPrefix(m_name, name);
  axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
  if(sidreGroup == nullptr)
  {
    return std::shared_ptr<Field>(nullptr);
  }

  std::string value;
  if(m_reader->getString(fullName, value))
  {
    sidreGroup->createViewString("value", value);
  }
  return addField(sidreGroup, axom::sidre::DataTypeId::CHAR8_STR_ID, fullName, name);
}

std::shared_ptr<Table> Table::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");

  if(m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format(
      "[Inlet] Table has already defined "
      "required value: {0}",
      m_sidreGroup->getName());

    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return shared_from_this();
  }

  if(isRequired)
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
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");

  if(!m_sidreGroup->hasView("required"))
  {
    return false;
  }
  axom::sidre::View* valueView = m_sidreGroup->getView("required");
  if(valueView == nullptr)
  {
    //TODO: is this possible after it says it has the view?
    return false;
  }
  int8 intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      " boolean value named {0}",
      m_sidreGroup->getName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return false;
  }

  return (bool)intValue;
}

std::shared_ptr<Table> Table::registerVerifier(std::function<bool()> lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Table "
                              "already set: {0}",
                              m_name));
  m_verifier = lambda;
  return shared_from_this();
}

bool Table::verify()
{
  bool verified = true;
  // Verify this Table
  if(m_verifier && !m_verifier())
  {
    verified = false;
    SLIC_WARNING(fmt::format("[Inlet] Table failed verification: {0}", m_name));
  }
  // Verify the child Fields of this Table
  for(auto field : m_fieldChildren)
  {
    if(!field.second->verify())
    {
      verified = false;
    }
  }
  // Verify the child Tables of this Table
  for(auto table : m_tableChildren)
  {
    if(!table.second->verify())
    {
      verified = false;
    }
  }

  return verified;
}

bool Table::hasChildTable(const std::string& tableName)
{
  return m_tableChildren.find(appendPrefix(m_name, tableName)) !=
    m_tableChildren.end();
}

bool Table::hasChildField(const std::string& fieldName)
{
  return m_fieldChildren.find(appendPrefix(m_name, fieldName)) !=
    m_fieldChildren.end();
}

bool Table::hasTable(const std::string& tableName)
{
  return static_cast<bool>(getTableInternal(tableName));
}

bool Table::hasField(const std::string& fieldName)
{
  return static_cast<bool>(getFieldInternal(fieldName));
}

std::shared_ptr<Table> Table::getTable(const std::string& tableName)
{
  auto table = getTableInternal(tableName);
  if(!table)
  {
    SLIC_WARNING(fmt::format("[Inlet] Table not found: {0}", tableName));
  }
  return table;
}

std::shared_ptr<Field> Table::getField(const std::string& fieldName)
{
  auto field = getFieldInternal(fieldName);
  if(!field)
  {
    SLIC_WARNING(fmt::format("[Inlet] Field not found: {0}", fieldName));
  }
  return field;
}

std::shared_ptr<Table> Table::getTableInternal(const std::string& tableName)
{
  std::string name = tableName;
  size_t found = name.find("/");
  auto currTable = shared_from_this();

  while(found != std::string::npos)
  {
    const std::string& currName = name.substr(0, found);
    if(currTable->hasChildTable(currName))
    {
      currTable =
        currTable->m_tableChildren[appendPrefix(currTable->m_name, currName)];
    }
    else
    {
      return nullptr;
    }
    name = name.substr(found + 1);
    found = name.find("/");
  }

  if(currTable->hasChildTable(name))
  {
    return currTable->m_tableChildren[appendPrefix(currTable->m_name, name)];
  }
  return nullptr;
}

std::shared_ptr<Field> Table::getFieldInternal(const std::string& fieldName)
{
  size_t found = fieldName.find_last_of("/");
  if(found == std::string::npos)
  {
    if(hasChildField(fieldName))
    {
      return m_fieldChildren[appendPrefix(m_name, fieldName)];
    }
  }
  else
  {
    const std::string& name = fieldName.substr(found + 1);
    auto table = getTableInternal(fieldName.substr(0, found));
    if(table && table->hasChildField(name))
    {
      return table->m_fieldChildren[appendPrefix(table->m_name, name)];
    }
  }
  return nullptr;
}

std::string Table::name() { return m_name; }

std::unordered_map<std::string, std::shared_ptr<Table>> Table::getChildTables()
{
  return m_tableChildren;
}

std::unordered_map<std::string, std::shared_ptr<Field>> Table::getChildFields()
{
  return m_fieldChildren;
}

}  // end namespace inlet
}  // end namespace axom
