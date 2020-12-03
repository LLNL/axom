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
const std::string Table::ARRAY_GROUP_NAME = "_inlet_array";
const std::string Table::ARRAY_INDICIES_VIEW_NAME = "_inlet_array_indices";

Table& Table::addTable(const std::string& name, const std::string& description)
{
  // Create intermediate Tables if they don't already exist
  std::string currName = name;
  size_t found = currName.find("/");
  auto currTable = this;

  while(found != std::string::npos)
  {
    const std::string currTableName =
      appendPrefix(currTable->m_name, currName.substr(0, found));
    // The current table will prepend its own prefix - just pass the basename
    if(!currTable->hasChild<Table>(currName.substr(0, found)))
    {
      // Will the copy always be elided here with a move ctor
      // or do we need std::piecewise_construct/std::forward_as_tuple?
      const auto& emplace_result = currTable->m_tableChildren.emplace(
        currTableName,
        cpp11_compat::make_unique<Table>(currTableName,
                                         "",
                                         m_reader,
                                         m_sidreRootGroup,
                                         m_docEnabled));
      // emplace_result is a pair whose first element is an iterator to the inserted element
      currTable = emplace_result.first->second.get();
    }
    else
    {
      currTable = currTable->m_tableChildren[currTableName].get();
    }
    currName = currName.substr(found + 1);
    found = currName.find("/");
  }

  const std::string currTableName =
    appendPrefix(currTable->m_name, currName.substr(0, found));

  if(!currTable->hasChild<Table>(currName))
  {
    const auto& emplace_result = currTable->m_tableChildren.emplace(
      currTableName,
      cpp11_compat::make_unique<Table>(currTableName,
                                       description,
                                       m_reader,
                                       m_sidreRootGroup,
                                       m_docEnabled));
    currTable = emplace_result.first->second.get();
  }
  else
  {
    currTable = currTable->m_tableChildren[currTableName].get();
  }

  return *currTable;
}

std::vector<std::pair<std::string, std::string>> Table::arrayIndicesWithPaths(
  const std::string& name) const
{
  std::vector<std::pair<std::string, std::string>> result;
  if(!m_sidreGroup->hasView(ARRAY_INDICIES_VIEW_NAME))
  {
    SLIC_ERROR(fmt::format(
      "[Inlet] Table '{0}' does not contain an array of user-defined objects",
      m_name));
  }
  const auto view = m_sidreGroup->getView(ARRAY_INDICIES_VIEW_NAME);
  const int* array = view->getArray();
  // Need to go up one level because this is an _inlet_array group
  const auto pos = m_name.find_last_of("/");
  const std::string baseName = m_name.substr(0, pos);
  for(int i = 0; i < view->getNumElements(); i++)
  {
    const auto indexLabel = std::to_string(array[i]);
    // The base name reflects the structure of the actual data
    // and is used for the reader call
    auto fullPath = appendPrefix(baseName, indexLabel);
    fullPath = appendPrefix(fullPath, name);
    // e.g. full_path could be foo/1/bar for field "bar" at index 1 of array "foo"
    result.push_back({indexLabel, fullPath});
  }
  return result;
}

Verifiable<Table>& Table::addBoolArray(const std::string& name,
                                       const std::string& description)
{
  return addPrimitiveArray<bool>(name, description);
}

Verifiable<Table>& Table::addIntArray(const std::string& name,
                                      const std::string& description)
{
  return addPrimitiveArray<int>(name, description);
}

Verifiable<Table>& Table::addDoubleArray(const std::string& name,
                                         const std::string& description)
{
  return addPrimitiveArray<double>(name, description);
}

Verifiable<Table>& Table::addStringArray(const std::string& name,
                                         const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description);
}

Table& Table::addGenericArray(const std::string& name,
                              const std::string& description)
{
  if(m_sidreGroup->hasView(ARRAY_INDICIES_VIEW_NAME))
  {
    SLIC_ERROR(
      fmt::format("[Inlet] Adding array of structs to array of structs {0} is "
                  "not supported",
                  m_name));
  }
  auto& table = addTable(appendPrefix(name, ARRAY_GROUP_NAME), description);
  std::vector<int> indices;
  const std::string& fullName = appendPrefix(m_name, name);
  if(m_reader.getArrayIndices(fullName, indices))
  {
    // This is how an array of user-defined type is differentiated
    // from an array of primitives - the tables have to be allocated
    // before they are populated as we don't know the schema of the
    // generic type yet
    auto view =
      table.m_sidreGroup->createViewAndAllocate(ARRAY_INDICIES_VIEW_NAME,
                                                axom::sidre::INT_ID,
                                                indices.size());
    int* raw_array = view->getArray();
    std::copy(indices.begin(), indices.end(), raw_array);
    for(const auto idx : indices)
    {
      table.addTable(std::to_string(idx), description);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("[Inlet] Array {0} not found.", fullName));
  }
  return table;
}

axom::sidre::Group* Table::createSidreGroup(const std::string& name,
                                            const std::string& description)
{
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

Field& Table::addField(axom::sidre::Group* sidreGroup,
                       axom::sidre::DataTypeId type,
                       const std::string& fullName,
                       const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currTable = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Tables (if not present) before adding the field
    currTable = &addTable(name.substr(0, found));
  }
  const auto& emplace_result = currTable->m_fieldChildren.emplace(
    fullName,
    cpp11_compat::make_unique<Field>(sidreGroup,
                                     m_sidreRootGroup,
                                     type,
                                     m_docEnabled));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

Function& Table::addFunctionInternal(axom::sidre::Group* sidreGroup,
                                     FunctionVariant&& func,
                                     const std::string& fullName,
                                     const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currTable = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Tables (if not present) before adding the field
    currTable = &addTable(name.substr(0, found));
  }
  const auto& emplace_result = currTable->m_functionChildren.emplace(
    fullName,
    cpp11_compat::make_unique<Function>(sidreGroup,
                                        m_sidreRootGroup,
                                        std::move(func),
                                        m_docEnabled));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

template <typename T, typename SFINAE>
VerifiableScalar& Table::addPrimitive(const std::string& name,
                                      const std::string& description,
                                      bool forArray,
                                      T val,
                                      const std::string& pathOverride)
{
  if(m_sidreGroup->hasView(ARRAY_INDICIES_VIEW_NAME))
  {
    // If it has indices, we're adding a primitive field to an array
    // of structs, so we need to iterate over the subtables
    // corresponding to elements of the array
    std::vector<std::reference_wrapper<VerifiableScalar>> fields;
    for(const auto& indexPath : arrayIndicesWithPaths(name))
    {
      // Add a primitive to an array element (which is a struct)
      fields.push_back(
        getTable(indexPath.first)
          .addPrimitive<T>(name, description, forArray, val, indexPath.second));
    }
    // Create an aggregate field so requirements can be collectively imposed
    // on all elements of the array
    m_aggregate_fields.emplace_back(std::move(fields));

    // Remove when C++17 is available
    return m_aggregate_fields.back();
  }
  else
  {
    // Otherwise actually add a Field
    std::string fullName = appendPrefix(m_name, name);
    axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
    SLIC_ERROR_IF(
      sidreGroup == nullptr,
      fmt::format("Failed to create Sidre group with name {0}", fullName));
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    auto typeId = addPrimitiveHelper(sidreGroup, lookupPath, forArray, val);
    return addField(sidreGroup, typeId, fullName, name);
  }
}

template <>
axom::sidre::DataTypeId Table::addPrimitiveHelper<bool>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  bool val)
{
  if(forArray || m_reader.getBool(lookupPath, val))
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
  if(forArray || m_reader.getInt(lookupPath, val))
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
  if(forArray || m_reader.getDouble(lookupPath, val))
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
  if(forArray || m_reader.getString(lookupPath, val))
  {
    sidreGroup->createViewString("value", val);
  }
  return axom::sidre::DataTypeId::CHAR8_STR_ID;
}

template <typename T, typename SFINAE>
Verifiable<Table>& Table::addPrimitiveArray(const std::string& name,
                                            const std::string& description,
                                            const std::string& pathOverride)
{
  if(m_sidreGroup->hasView(ARRAY_INDICIES_VIEW_NAME))
  {
    // Adding an array of primitive field to an array of structs
    std::vector<std::reference_wrapper<Verifiable>> tables;
    // Iterate over each element and forward the call to addPrimitiveArray
    for(const auto& indexPath : arrayIndicesWithPaths(name))
    {
      tables.push_back(
        getTable(indexPath.first)
          .addPrimitiveArray<T>(name, description, indexPath.second));
    }

    m_aggregate_tables.emplace_back(std::move(tables));

    // Remove when C++17 is available
    return m_aggregate_tables.back();
  }
  else
  {
    // "base case", create a table for the field and fill it in with the helper
    auto& table = addTable(appendPrefix(name, ARRAY_GROUP_NAME), description);
    const std::string& fullName = appendPrefix(m_name, name);
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    addPrimitiveArrayHelper<T>(table, lookupPath);
    return table;
  }
}

template <>
void Table::addPrimitiveArrayHelper<bool>(Table& table,
                                          const std::string& lookupPath)
{
  std::unordered_map<int, bool> map;
  if(m_reader.getBoolMap(lookupPath, map))
  {
    for(const auto& p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("[Inlet] Bool array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<int>(Table& table,
                                         const std::string& lookupPath)
{
  std::unordered_map<int, int> map;
  if(m_reader.getIntMap(lookupPath, map))
  {
    for(const auto& p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("[Inlet] Int array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<double>(Table& table,
                                            const std::string& lookupPath)
{
  std::unordered_map<int, double> map;
  if(m_reader.getDoubleMap(lookupPath, map))
  {
    for(const auto& p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("[Inlet] Double array {0} not found.", lookupPath));
  }
}

template <>
void Table::addPrimitiveArrayHelper<std::string>(Table& table,
                                                 const std::string& lookupPath)
{
  std::unordered_map<int, std::string> map;
  if(m_reader.getStringMap(lookupPath, map))
  {
    for(const auto& p : map)
    {
      table.addPrimitive(std::to_string(p.first), "", true, p.second);
    }
  }
  else
  {
    SLIC_WARNING(fmt::format("[Inlet] String array {0} not found.", lookupPath));
  }
}

Verifiable<Function>& Table::addFunction(const std::string& name,
                                         const FunctionType ret_type,
                                         const std::vector<FunctionType>& arg_types,
                                         const std::string& description,
                                         const std::string& pathOverride)
{
  if(m_sidreGroup->hasView(ARRAY_INDICIES_VIEW_NAME))
  {
    // If it has indices, we're adding a primitive field to an array
    // of structs, so we need to iterate over the subtables
    // corresponding to elements of the array
    std::vector<std::reference_wrapper<Verifiable<Function>>> funcs;
    for(const auto& indexPath : arrayIndicesWithPaths(name))
    {
      // Add a primitive to an array element (which is a struct)
      funcs.push_back(
        getTable(indexPath.first)
          .addFunction(name, ret_type, arg_types, description, indexPath.second));
    }
    // Create an aggregate field so requirements can be collectively imposed
    // on all elements of the array
    m_aggregate_funcs.emplace_back(std::move(funcs));

    // Remove when C++17 is available
    return m_aggregate_funcs.back();
  }
  else
  {
    // Otherwise actually add a Field
    std::string fullName = appendPrefix(m_name, name);
    axom::sidre::Group* sidreGroup = createSidreGroup(fullName, description);
    SLIC_ERROR_IF(
      sidreGroup == nullptr,
      fmt::format("Failed to create Sidre group with name {0}", fullName));
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    auto func = m_reader.getFunction(lookupPath, ret_type, arg_types);
    return addFunctionInternal(sidreGroup, std::move(func), fullName, name);
  }
}

Proxy Table::operator[](const std::string& name) const
{
  const bool has_table = hasTable(name);
  const bool has_field = hasField(name);
  const bool has_func = hasFunction(name);

  // Ambiguous case - both a table and field exist with the same name
  if((has_table && has_field) || (has_field && has_func) ||
     (has_table && has_func))
  {
    const std::string msg = fmt::format(
      "[Inlet] Ambiguous lookup - more than one of a table/field/function with "
      "name {0} exist",
      name);
    SLIC_ERROR(msg);
    return Proxy();
  }

  else if(has_table)
  {
    return Proxy(getTable(name));
  }

  else if(has_field)
  {
    return Proxy(getField(name));
  }

  else if(has_func)
  {
    return Proxy(getFunction(name));
  }

  // Neither exists
  else
  {
    std::string msg =
      fmt::format("[Inlet] No table, field, or function with name {0} exists",
                  name);
    SLIC_ERROR(msg);
    return Proxy();
  }
}

Table& Table::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");
  setRequired(*m_sidreGroup, *m_sidreRootGroup, isRequired);
  return *this;
}

bool Table::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");
  return checkIfRequired(*m_sidreGroup, *m_sidreRootGroup);
}

Table& Table::registerVerifier(std::function<bool(const Table&)> lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Table "
                              "already set: {0}",
                              m_name));
  m_verifier = lambda;
  return *this;
}

bool Table::verify() const
{
  bool verified = true;
  // If this table was required, make sure something was defined in it
  verified &=
    verifyRequired(*m_sidreGroup, m_sidreGroup->getNumGroups() > 0, "Table");
  // Verify this Table if a lambda was configured
  if(m_verifier && !m_verifier(*this))
  {
    verified = false;
    SLIC_WARNING(fmt::format("[Inlet] Table failed verification: {0}", m_name));
  }
  // Verify the child Fields of this Table
  for(const auto& field : m_fieldChildren)
  {
    if(!field.second->verify())
    {
      verified = false;
    }
  }
  // Verify the child Tables of this Table
  for(const auto& table : m_tableChildren)
  {
    if(!table.second->verify())
    {
      verified = false;
    }
  }

  // Verify the child Functions of this Table
  for(const auto& function : m_functionChildren)
  {
    if(!function.second->verify())
    {
      verified = false;
    }
  }

  return verified;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Table>> Table::*Table::getChildren()
{
  return &Table::m_tableChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Field>> Table::*Table::getChildren()
{
  return &Table::m_fieldChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Function>> Table::*Table::getChildren()
{
  return &Table::m_functionChildren;
}

template <typename T>
bool Table::hasChild(const std::string& childName) const
{
  const auto& children = this->*getChildren<T>();
  return children.find(appendPrefix(m_name, childName)) != children.end();
}

template <typename T>
T* Table::getChildInternal(const std::string& childName) const
{
  std::string name = childName;
  size_t found = name.find("/");
  auto currTable = this;

  while(found != std::string::npos)
  {
    const std::string& currName = name.substr(0, found);
    if(currTable->hasChild<Table>(currName))
    {
      currTable =
        currTable->m_tableChildren.at(appendPrefix(currTable->m_name, currName))
          .get();
    }
    else
    {
      return nullptr;
    }
    name = name.substr(found + 1);
    found = name.find("/");
  }

  if(currTable->hasChild<T>(name))
  {
    const auto& children = currTable->*getChildren<T>();
    return children.at(appendPrefix(currTable->m_name, name)).get();
  }
  return nullptr;
}

bool Table::hasTable(const std::string& tableName) const
{
  return static_cast<bool>(getChildInternal<Table>(tableName));
}

bool Table::hasField(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Field>(fieldName));
}

bool Table::hasFunction(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Function>(fieldName));
}

Table& Table::getTable(const std::string& tableName) const
{
  auto table = getChildInternal<Table>(tableName);
  if(!table)
  {
    SLIC_ERROR(fmt::format("[Inlet] Table not found: {0}", tableName));
  }
  return *table;
}

Field& Table::getField(const std::string& fieldName) const
{
  auto field = getChildInternal<Field>(fieldName);
  if(!field)
  {
    SLIC_ERROR(fmt::format("[Inlet] Field not found: {0}", fieldName));
  }
  return *field;
}

Function& Table::getFunction(const std::string& funcName) const
{
  auto func = getChildInternal<Function>(funcName);
  if(!func)
  {
    SLIC_ERROR(fmt::format("[Inlet] Function not found: {0}", funcName));
  }
  return *func;
}

std::string Table::name() const { return m_name; }

bool Table::contains(const std::string& name) const
{
  if(auto table = getChildInternal<Table>(name))
  {
    // call operator bool on the table itself
    return static_cast<bool>(*table);
  }
  else if(auto field = getChildInternal<Field>(name))
  {
    // call operator bool on the field itself
    return static_cast<bool>(*field);
  }
  else if(auto function = getChildInternal<Function>(name))
  {
    // call operator bool on the function itself
    return static_cast<bool>(*function);
  }
  return false;
}

Table::operator bool() const
{
  // If there are no child tables, it must have child fields
  if(m_tableChildren.empty())
  {
    return !m_fieldChildren.empty();
  }

  // Otherwise we have to recurse and check child tables
  return std::any_of(m_tableChildren.begin(),
                     m_tableChildren.end(),
                     [](const decltype(m_tableChildren)::value_type& entry) {
                       return static_cast<bool>(entry.second);
                     });
}

const std::unordered_map<std::string, std::unique_ptr<Table>>&
Table::getChildTables() const
{
  return m_tableChildren;
}

const std::unordered_map<std::string, std::unique_ptr<Field>>&
Table::getChildFields() const
{
  return m_fieldChildren;
}

}  // end namespace inlet
}  // end namespace axom
