// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Table.hpp"

#include "axom/slic.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/Proxy.hpp"

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
    // The current table will prepend its own prefix - just pass the basename
    if(!currTable->hasChildTable(currName.substr(0, found)))
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

std::vector<std::pair<std::string, std::string>> Table::arrayIndicesWithPaths(
  const std::string& name)
{
  std::vector<std::pair<std::string, std::string>> result;
  if(!m_sidreGroup->hasView("_inlet_array_indices"))
  {
    SLIC_ERROR(fmt::format(
      "Table '{0}' does not contain an array of user-defined objects",
      m_name));
  }
  auto view = m_sidreGroup->getView("_inlet_array_indices");
  int* array = view->getArray();
  // Need to go up one level because this is an _inlet_array group
  auto pos = m_name.find_last_of("/");
  std::string baseName = m_name.substr(0, pos);
  for(int i = 0; i < view->getNumElements(); i++)
  {
    auto indexLabel = std::to_string(array[i]);
    // The base name reflects the structure of the actual data
    // and is used for the reader call
    auto fullPath = appendPrefix(baseName, indexLabel);
    fullPath = appendPrefix(fullPath, name);
    // e.g. full_path could be foo/1/bar for field "bar" at index 1 of array "foo"
    result.push_back({indexLabel, fullPath});
  }
  return result;
}

std::shared_ptr<Verifiable> Table::addBoolArray(const std::string& name,
                                                const std::string& description)
{
  return addPrimitiveArray<bool>(name, description);
}

std::shared_ptr<Verifiable> Table::addIntArray(const std::string& name,
                                               const std::string& description)
{
  return addPrimitiveArray<int>(name, description);
}

std::shared_ptr<Verifiable> Table::addDoubleArray(const std::string& name,
                                                  const std::string& description)
{
  return addPrimitiveArray<double>(name, description);
}

std::shared_ptr<Verifiable> Table::addStringArray(const std::string& name,
                                                  const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description);
}

std::shared_ptr<Table> Table::addGenericArray(const std::string& name,
                                              const std::string& description)
{
  auto table = addTable(appendPrefix(name, "_inlet_array"), description);
  std::vector<int> indices;
  const std::string& fullName = appendPrefix(m_name, name);
  if(m_reader->getArrayIndices(fullName, indices))
  {
    // This is how an array of user-defined type is differentiated
    // from an array of primitives - the tables have to be allocated
    // before they are populated as we don't know the schema of the
    // generic type yet
    auto view =
      table->m_sidreGroup->createViewAndAllocate("_inlet_array_indices",
                                                 axom::sidre::INT_ID,
                                                 indices.size());
    int* raw_array = view->getArray();
    std::copy(indices.begin(), indices.end(), raw_array);
    for(const auto idx : indices)
    {
      table->addTable(std::to_string(idx), description);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("Array {0} not found.", fullName));
  }
  return table;
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

template <typename T, typename SFINAE>
std::shared_ptr<VerifiableScalar> Table::addPrimitive(
  const std::string& name,
  const std::string& description,
  bool forArray,
  T val,
  const std::string& pathOverride)
{
  if(m_sidreGroup->hasView("_inlet_array_indices"))
  {
    std::vector<std::shared_ptr<VerifiableScalar>> fields;
    for(const auto& indexPath : arrayIndicesWithPaths(name))
    {
      fields.push_back(
        getTable(indexPath.first)
          ->addPrimitive<T>(name, description, forArray, val, indexPath.second));
    }
    return std::make_shared<AggregateField>(std::move(fields));
  }
  else
  {
    std::string fullName = appendPrefix(m_name, name);
    axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
    if(sidreGroup == nullptr)
    {
      return std::shared_ptr<Field>(nullptr);
    }

    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    auto typeId = addPrimitiveHelper(sidreGroup, lookupPath, forArray, val);
    return addField(sidreGroup, typeId, fullName, name);
  }
}

// Explicit instantiations
template std::shared_ptr<VerifiableScalar> Table::addPrimitive<bool>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  bool val,
  const std::string& pathOverride);

template std::shared_ptr<VerifiableScalar> Table::addPrimitive<int>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  int val,
  const std::string& pathOverride);

template std::shared_ptr<VerifiableScalar> Table::addPrimitive<double>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  double val,
  const std::string& pathOverride);

template std::shared_ptr<VerifiableScalar> Table::addPrimitive<std::string>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  std::string val,
  const std::string& pathOverride);

template <>
axom::sidre::DataTypeId Table::addPrimitiveHelper<bool>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  bool val)
{
  if(forArray || m_reader->getBool(lookupPath, val))
  {
    sidreGroup->createViewScalar("value", val ? int8(1) : int8(0));
  }
  return axom::sidre::DataTypeId::INT8_ID;
}

template <>
axom::sidre::DataTypeId Table::addPrimitiveHelper<int>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  int val)
{
  if(forArray || m_reader->getInt(lookupPath, val))
  {
    sidreGroup->createViewScalar("value", val);
  }
  return axom::sidre::DataTypeId::INT_ID;
}

template <>
axom::sidre::DataTypeId Table::addPrimitiveHelper<double>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  double val)
{
  if(forArray || m_reader->getDouble(lookupPath, val))
  {
    sidreGroup->createViewScalar("value", val);
  }
  return axom::sidre::DataTypeId::DOUBLE_ID;
}

template <>
axom::sidre::DataTypeId Table::addPrimitiveHelper<std::string>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  std::string val)
{
  if(forArray || m_reader->getString(lookupPath, val))
  {
    sidreGroup->createViewString("value", val);
  }
  return axom::sidre::DataTypeId::CHAR8_STR_ID;
}

template <typename T, typename SFINAE>
std::shared_ptr<Verifiable> Table::addPrimitiveArray(
  const std::string& name,
  const std::string& description,
  const std::string& pathOverride)
{
  if(m_sidreGroup->hasView("_inlet_array_indices"))
  {
    std::vector<std::shared_ptr<Verifiable>> tables;
    for(const auto& indexPath : arrayIndicesWithPaths(name))
    {
      tables.push_back(
        getTable(indexPath.first)
          ->addPrimitiveArray<T>(name, description, indexPath.second));
    }
    return std::make_shared<AggregateTable>(std::move(tables));
  }
  else
  {
    auto table = addTable(appendPrefix(name, "_inlet_array"), description);
    const std::string& fullName = appendPrefix(m_name, name);
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    addPrimitiveArrayHelper<T>(*table, lookupPath);
    return table;
  }
}

// Explicit instantiations
template std::shared_ptr<Verifiable> Table::addPrimitiveArray<bool>(
  const std::string& name,
  const std::string& description,
  const std::string& pathOverride);

template std::shared_ptr<Verifiable> Table::addPrimitiveArray<int>(
  const std::string& name,
  const std::string& description,
  const std::string& pathOverride);

template std::shared_ptr<Verifiable> Table::addPrimitiveArray<double>(
  const std::string& name,
  const std::string& description,
  const std::string& pathOverride);

template std::shared_ptr<Verifiable> Table::addPrimitiveArray<std::string>(
  const std::string& name,
  const std::string& description,
  const std::string& pathOverride);

template <>
void Table::addPrimitiveArrayHelper<bool>(Table& table,
                                          const std::string& lookupPath)
{
  std::unordered_map<int, bool> map;
  if(m_reader->getBoolMap(lookupPath, map))
  {
    for(auto p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("Bool array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<int>(Table& table,
                                         const std::string& lookupPath)
{
  std::unordered_map<int, int> map;
  if(m_reader->getIntMap(lookupPath, map))
  {
    for(auto p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("Int array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<double>(Table& table,
                                            const std::string& lookupPath)
{
  std::unordered_map<int, double> map;
  if(m_reader->getDoubleMap(lookupPath, map))
  {
    for(auto p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("Double array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<std::string>(Table& table,
                                                 const std::string& lookupPath)
{
  std::unordered_map<int, std::string> map;
  if(m_reader->getStringMap(lookupPath, map))
  {
    for(auto p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("String array {0} not found.", lookupPath));
  }
}

Proxy Table::operator[](const std::string& name)
{
  auto has_table = hasTable(name);
  auto has_field = hasField(name);

  // Ambiguous case - both a table and field exist with the same name
  if(has_table && has_field)
  {
    std::string msg = fmt::format(
      "[Inlet] Ambiguous lookup - both a table and field with name {0} exist",
      name);
    SLIC_ERROR(msg);
    return Proxy();
  }

  else if(has_table)
  {
    return Proxy(*getTable(name));
  }

  else if(has_field)
  {
    return Proxy(*getField(name));
  }

  // Neither exists
  else
  {
    std::string msg =
      fmt::format("[Inlet] Neither a table nor a field with name {0} exist",
                  name);
    SLIC_ERROR(msg);
    return Proxy();
  }
}

std::shared_ptr<Verifiable> Table::required(bool isRequired)
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

bool Table::isRequired()
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

std::shared_ptr<Verifiable> Table::registerVerifier(std::function<bool(Table&)> lambda)
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
  // If this table was required, make sure soemething was defined in it
  if(m_sidreGroup->hasView("required"))
  {
    int8 required = m_sidreGroup->getView("required")->getData();
    if(required && m_sidreGroup->getNumGroups() == 0)
    {
      std::string msg = fmt::format(
        "[Inlet] Required Table not "
        "specified: {0}",
        m_sidreGroup->getPathName());
      SLIC_WARNING(msg);
      verified = false;
    }
  }
  // Verify this Table if a lambda was configured
  if(m_verifier && !m_verifier(*this))
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

bool AggregateTable::verify()
{
  return std::all_of(
    m_tables.begin(),
    m_tables.end(),
    [](std::shared_ptr<Verifiable>& table) { return table->verify(); });
}

std::shared_ptr<Verifiable> AggregateTable::required(bool isRequired)
{
  for(auto& table : m_tables)
  {
    table->required(isRequired);
  }
  return shared_from_this();
}

bool AggregateTable::isRequired()
{
  for(auto& table : m_tables)
  {
    if(table->isRequired())
    {
      return true;
    }
  }
  return false;
}

std::shared_ptr<Verifiable> AggregateTable::registerVerifier(
  std::function<bool(Table&)> lambda)
{
  for(auto& table : m_tables)
  {
    table->registerVerifier(lambda);
  }
  return shared_from_this();
}

}  // end namespace inlet
}  // end namespace axom
