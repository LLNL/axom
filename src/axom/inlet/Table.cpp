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
const std::string Table::CONTAINER_GROUP_NAME = "_inlet_container";
const std::string Table::CONTAINER_INDICES_NAME = "_inlet_container_indices";

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
    if(!currTable->hasChildTable(currName.substr(0, found)))
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

  if(!currTable->hasChildTable(currName))
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

std::vector<std::string> Table::containerIndices() const
{
  std::vector<std::string> indices;
  // If it's in the view, they are integer indices of an array
  if(m_sidreGroup->hasView(CONTAINER_INDICES_NAME))
  {
    const auto view = m_sidreGroup->getView(CONTAINER_INDICES_NAME);
    const int* array = view->getArray();
    const auto num_elements = view->getNumElements();
    indices.reserve(num_elements);
    for(int i = 0; i < num_elements; i++)
    {
      indices.push_back(std::to_string(array[i]));
    }
  }
  // Otherwise, they're string indices of a dict
  // need to iterate over group
  else if(m_sidreGroup->hasGroup(CONTAINER_INDICES_NAME))
  {
    auto group = m_sidreGroup->getGroup(CONTAINER_INDICES_NAME);
    indices.reserve(group->getNumViews());
    for(auto idx = group->getFirstValidViewIndex(); sidre::indexIsValid(idx);
        idx = group->getNextValidViewIndex(idx))
    {
      indices.push_back(group->getView(idx)->getString());
    }
  }
  else
  {
    SLIC_ERROR(
      fmt::format("[Inlet] Table '{0}' does not contain an array or dict of "
                  "user-defined objects",
                  m_name));
  }
  return indices;
}

std::vector<std::pair<std::string, std::string>> Table::containerIndicesWithPaths(
  const std::string& name) const
{
  std::vector<std::pair<std::string, std::string>> result;
  // Need to go up one level because this is an _inlet_array group
  const auto pos = m_name.find_last_of("/");
  const std::string baseName = m_name.substr(0, pos);
  for(const auto& indexLabel : containerIndices())
  {
    // The base name reflects the structure of the actual data
    // and is used for the reader call
    auto fullPath = appendPrefix(baseName, indexLabel);
    fullPath = appendPrefix(fullPath, name);
    // e.g. full_path could be foo/1/bar for field "bar" at index 1 of array "foo"
    result.push_back({indexLabel, fullPath});
  }
  return result;
}

Verifiable& Table::addBoolArray(const std::string& name,
                                const std::string& description)
{
  return addPrimitiveArray<bool>(name, description);
}

Verifiable& Table::addIntArray(const std::string& name,
                               const std::string& description)
{
  return addPrimitiveArray<int>(name, description);
}

Verifiable& Table::addDoubleArray(const std::string& name,
                                  const std::string& description)
{
  return addPrimitiveArray<double>(name, description);
}

Verifiable& Table::addStringArray(const std::string& name,
                                  const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description);
}

Table& Table::addGenericArray(const std::string& name,
                              const std::string& description)
{
  if(isGenericContainer())
  {
    SLIC_ERROR(fmt::format(
      "[Inlet] Adding array of container to array of structs {0} is "
      "not supported",
      m_name));
  }
  auto& table = addTable(appendPrefix(name, CONTAINER_GROUP_NAME), description);
  std::vector<int> indices;
  const std::string& fullName = appendPrefix(m_name, name);
  if(m_reader.getIndices(fullName, indices))
  {
    // This is how an array of user-defined type is differentiated
    // from an array of primitives - the tables have to be allocated
    // before they are populated as we don't know the schema of the
    // generic type yet
    auto view = table.m_sidreGroup->createViewAndAllocate(CONTAINER_INDICES_NAME,
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

Verifiable& Table::addBoolDict(const std::string& name,
                               const std::string& description)
{
  return addPrimitiveArray<bool>(name, description, true);
}

Verifiable& Table::addIntDict(const std::string& name,
                              const std::string& description)
{
  return addPrimitiveArray<int>(name, description, true);
}

Verifiable& Table::addDoubleDict(const std::string& name,
                                 const std::string& description)
{
  return addPrimitiveArray<double>(name, description, true);
}

Verifiable& Table::addStringDict(const std::string& name,
                                 const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description, true);
}

Table& Table::addGenericDict(const std::string& name,
                             const std::string& description)
{
  if(isGenericContainer())
  {
    SLIC_ERROR(fmt::format(
      "[Inlet] Adding dict of structs to container of structs {0} is "
      "not supported",
      m_name));
  }
  auto& table = addTable(appendPrefix(name, CONTAINER_GROUP_NAME), description);
  std::vector<std::string> indices;
  const std::string& fullName = appendPrefix(m_name, name);
  if(m_reader.getIndices(fullName, indices))
  {
    auto group = table.m_sidreGroup->createGroup(CONTAINER_INDICES_NAME,
                                                 /* list_format = */ true);
    for(const auto idx : indices)
    {
      table.addTable(idx, description);
      group->createViewString("", idx);
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

template <typename T, typename SFINAE>
VerifiableScalar& Table::addPrimitive(const std::string& name,
                                      const std::string& description,
                                      bool forArray,
                                      T val,
                                      const std::string& pathOverride)
{
  if(isGenericContainer())
  {
    // If it has indices, we're adding a primitive field to an array
    // of structs, so we need to iterate over the subtables
    // corresponding to elements of the array
    std::vector<std::reference_wrapper<VerifiableScalar>> fields;
    for(const auto& indexPath : containerIndicesWithPaths(name))
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

// Explicit instantiations
template VerifiableScalar& Table::addPrimitive<bool>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  bool val,
  const std::string& pathOverride);

template VerifiableScalar& Table::addPrimitive<int>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  int val,
  const std::string& pathOverride);

template VerifiableScalar& Table::addPrimitive<double>(
  const std::string& name,
  const std::string& description,
  bool forArray,
  double val,
  const std::string& pathOverride);

template VerifiableScalar& Table::addPrimitive<std::string>(
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

namespace detail
{
/*!
  *****************************************************************************
  * \brief Adds the contents of an array to the table
  *****************************************************************************
  */
template <typename T>
void registerContainer(Table& table, const std::unordered_map<int, T>& container)
{
  for(const auto& entry : container)
  {
    table.addPrimitive(std::to_string(entry.first), "", true, entry.second);
  }
}

/*!
  *****************************************************************************
  * \brief Adds the contents of a dict to the table
  *****************************************************************************
  */
template <typename T>
void registerContainer(Table& table,
                       const std::unordered_map<std::string, T>& container)
{
  for(const auto& entry : container)
  {
    table.addPrimitive(entry.first, "", true, entry.second);
  }
}
/*!
 *****************************************************************************
 * \brief Implementation helper for adding primitive arrays
 * 
 * \note Structs are used for partial template specializations
 *****************************************************************************
 */
template <typename Key, typename Primitive>
struct PrimitiveArrayHelper
{ };

template <typename Key>
struct PrimitiveArrayHelper<Key, bool>
{
  /*!
   *****************************************************************************
   * \brief Finalizes the creation of a container
   * \param [inout] table The table to add the container to
   * \param [in] reader The Reader object to read the container from
   * \param [in] lookupPath The path within the input file to the container
   *****************************************************************************
   */
  PrimitiveArrayHelper(Table& table, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, bool> map;
    if(!reader.getBoolMap(lookupPath, map))
    {
      SLIC_WARNING(
        fmt::format("[Inlet] Bool container {0} not found.", lookupPath));
    }
    registerContainer(table, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, int>
{
  PrimitiveArrayHelper(Table& table, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, int> map;
    if(!reader.getIntMap(lookupPath, map))
    {
      SLIC_WARNING(
        fmt::format("[Inlet] Int container {0} not found.", lookupPath));
    }
    registerContainer(table, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, double>
{
  PrimitiveArrayHelper(Table& table, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, double> map;
    if(!reader.getDoubleMap(lookupPath, map))
    {
      SLIC_WARNING(
        fmt::format("[Inlet] Double container {0} not found.", lookupPath));
    }
    registerContainer(table, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, std::string>
{
  PrimitiveArrayHelper(Table& table, Reader& reader, const std::string& lookupPath)
  {
    std::unordered_map<Key, std::string> map;
    if(!reader.getStringMap(lookupPath, map))
    {
      SLIC_WARNING(
        fmt::format("[Inlet] String container {0} not found.", lookupPath));
    }
    registerContainer(table, map);
  }
};

}  // end namespace detail

template <typename T, typename SFINAE>
Verifiable& Table::addPrimitiveArray(const std::string& name,
                                     const std::string& description,
                                     const bool isDict,
                                     const std::string& pathOverride)
{
  if(isGenericContainer())
  {
    // Adding an array of primitive field to an array of structs
    std::vector<std::reference_wrapper<Verifiable>> tables;
    // Iterate over each element and forward the call to addPrimitiveArray
    for(const auto& indexPath : containerIndicesWithPaths(name))
    {
      tables.push_back(
        getTable(indexPath.first)
          .addPrimitiveArray<T>(name, description, isDict, indexPath.second));
    }

    m_aggregate_tables.emplace_back(std::move(tables));

    // Remove when C++17 is available
    return m_aggregate_tables.back();
  }
  else
  {
    // "base case", create a table for the field and fill it in with the helper
    auto& table = addTable(appendPrefix(name, CONTAINER_GROUP_NAME), description);
    const std::string& fullName = appendPrefix(m_name, name);
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    if(isDict)
    {
      detail::PrimitiveArrayHelper<std::string, T>(table, m_reader, lookupPath);
    }
    else
    {
      detail::PrimitiveArrayHelper<int, T>(table, m_reader, lookupPath);
    }
    return table;
  }
}

// Explicit instantiations
template Verifiable& Table::addPrimitiveArray<bool>(
  const std::string& name,
  const std::string& description,
  const bool isDict,
  const std::string& pathOverride);

template Verifiable& Table::addPrimitiveArray<int>(const std::string& name,
                                                   const std::string& description,
                                                   const bool isDict,
                                                   const std::string& pathOverride);

template Verifiable& Table::addPrimitiveArray<double>(
  const std::string& name,
  const std::string& description,
  const bool isDict,
  const std::string& pathOverride);

template Verifiable& Table::addPrimitiveArray<std::string>(
  const std::string& name,
  const std::string& description,
  const bool isDict,
  const std::string& pathOverride);

Proxy Table::operator[](const std::string& name) const
{
  const bool has_table = hasTable(name);
  const bool has_field = hasField(name);

  // Ambiguous case - both a table and field exist with the same name
  if(has_table && has_field)
  {
    const std::string msg = fmt::format(
      "[Inlet] Ambiguous lookup - both a table and field with name {0} exist",
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

Table& Table::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");

  if(m_sidreGroup->hasView("required"))
  {
    const std::string msg = fmt::format(
      "[Inlet] Table has already defined "
      "required value: {0}",
      m_sidreGroup->getName());

    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return *this;
  }

  if(isRequired)
  {
    m_sidreGroup->createViewScalar("required", (int8)1);
  }
  else
  {
    m_sidreGroup->createViewScalar("required", (int8)0);
  }

  return *this;
}

bool Table::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Table specific Sidre Datastore Group not set");

  if(!m_sidreGroup->hasView("required"))
  {
    return false;
  }
  const axom::sidre::View* valueView = m_sidreGroup->getView("required");
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
      m_sidreGroup->getName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return false;
  }

  return (bool)intValue;
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
  // If this table was required, make sure soemething was defined in it
  if(m_sidreGroup->hasView("required"))
  {
    int8 required = m_sidreGroup->getView("required")->getData();
    if(required && m_sidreGroup->getNumGroups() == 0)
    {
      const std::string msg = fmt::format(
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

  return verified;
}

bool Table::hasChildTable(const std::string& tableName) const
{
  return m_tableChildren.find(appendPrefix(m_name, tableName)) !=
    m_tableChildren.end();
}

bool Table::hasChildField(const std::string& fieldName) const
{
  return m_fieldChildren.find(appendPrefix(m_name, fieldName)) !=
    m_fieldChildren.end();
}

bool Table::hasTable(const std::string& tableName) const
{
  return static_cast<bool>(getTableInternal(tableName));
}

bool Table::hasField(const std::string& fieldName) const
{
  return static_cast<bool>(getFieldInternal(fieldName));
}

Table& Table::getTable(const std::string& tableName) const
{
  auto table = getTableInternal(tableName);
  if(!table)
  {
    SLIC_ERROR(fmt::format("[Inlet] Table not found: {0}", tableName));
  }
  return *table;
}

Field& Table::getField(const std::string& fieldName) const
{
  auto field = getFieldInternal(fieldName);
  if(!field)
  {
    SLIC_ERROR(fmt::format("[Inlet] Field not found: {0}", fieldName));
  }
  return *field;
}

Table* Table::getTableInternal(const std::string& tableName) const
{
  std::string name = tableName;
  size_t found = name.find("/");
  auto currTable = this;

  while(found != std::string::npos)
  {
    const std::string& currName = name.substr(0, found);
    if(currTable->hasChildTable(currName))
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

  if(currTable->hasChildTable(name))
  {
    return currTable->m_tableChildren.at(appendPrefix(currTable->m_name, name)).get();
  }
  return nullptr;
}

Field* Table::getFieldInternal(const std::string& fieldName) const
{
  const size_t found = fieldName.find_last_of("/");
  if(found == std::string::npos)
  {
    if(hasChildField(fieldName))
    {
      return m_fieldChildren.at(appendPrefix(m_name, fieldName)).get();
    }
  }
  else
  {
    const std::string& name = fieldName.substr(found + 1);
    auto table = getTableInternal(fieldName.substr(0, found));
    if(table && table->hasChildField(name))
    {
      return table->m_fieldChildren.at(appendPrefix(table->m_name, name)).get();
    }
  }
  return nullptr;
}

std::string Table::name() const { return m_name; }

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

bool AggregateTable::verify() const
{
  return std::all_of(m_tables.begin(),
                     m_tables.end(),
                     [](const Verifiable& table) { return table.verify(); });
}

AggregateTable& AggregateTable::required(bool isRequired)
{
  for(auto& table : m_tables)
  {
    table.get().required(isRequired);
  }
  return *this;
}

bool AggregateTable::isRequired() const
{
  return std::any_of(m_tables.begin(),
                     m_tables.end(),
                     [](const Verifiable& table) { return table.isRequired(); });
}

AggregateTable& AggregateTable::registerVerifier(
  std::function<bool(const Table&)> lambda)
{
  for(auto& table : m_tables)
  {
    table.get().registerVerifier(lambda);
  }
  return *this;
}

}  // end namespace inlet
}  // end namespace axom
