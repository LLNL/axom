// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Container.hpp"

#include "axom/slic.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/Proxy.hpp"

namespace axom
{
namespace inlet
{
template <typename Func>
void Container::forEachCollectionElement(Func&& func) const
{
  for(const auto& index : collectionIndices())
  {
    func(getContainer(detail::indexToString(index)));
  }
}

Container& Container::addContainer(const std::string& name,
                                   const std::string& description)
{
  // Create intermediate Containers if they don't already exist
  std::string currName = name;
  size_t found = currName.find("/");
  auto currContainer = this;

  while(found != std::string::npos)
  {
    const std::string currContainerName =
      appendPrefix(currContainer->m_name, currName.substr(0, found));
    // The current container will prepend its own prefix - just pass the basename
    if(!currContainer->hasChild<Container>(currName.substr(0, found)))
    {
      // Will the copy always be elided here with a move ctor
      // or do we need std::piecewise_construct/std::forward_as_tuple?
      const auto& emplace_result = currContainer->m_containerChildren.emplace(
        currContainerName,
        cpp11_compat::make_unique<Container>(currContainerName,
                                             "",
                                             m_reader,
                                             m_sidreRootGroup,
                                             m_unexpectedNames,
                                             m_docEnabled));
      // emplace_result is a pair whose first element is an iterator to the inserted element
      currContainer = emplace_result.first->second.get();
    }
    else
    {
      currContainer = currContainer->m_containerChildren[currContainerName].get();
    }
    currName = currName.substr(found + 1);
    found = currName.find("/");
  }

  const std::string currContainerName =
    appendPrefix(currContainer->m_name, currName.substr(0, found));

  if(!currContainer->hasChild<Container>(currName))
  {
    const auto& emplace_result = currContainer->m_containerChildren.emplace(
      currContainerName,
      cpp11_compat::make_unique<Container>(currContainerName,
                                           description,
                                           m_reader,
                                           m_sidreRootGroup,
                                           m_unexpectedNames,
                                           m_docEnabled));
    currContainer = emplace_result.first->second.get();
  }
  else
  {
    currContainer = currContainer->m_containerChildren[currContainerName].get();
  }

  return *currContainer;
}

Container& Container::addStruct(const std::string& name,
                                const std::string& description)
{
  auto& base_container = addContainer(name, description);
  for(Container& sub_container : m_nested_aggregates)
  {
    base_container.m_nested_aggregates.push_back(
      sub_container.addStruct(name, description));
  }
  if(isStructCollection())
  {
    for(const auto& index : collectionIndices())
    {
      base_container.m_nested_aggregates.push_back(
        getContainer(detail::indexToString(index)).addStruct(name, description));
    }
  }
  return base_container;
}

std::vector<VariantKey> Container::collectionIndices(bool trimAbsolute) const
{
  std::vector<VariantKey> indices;
  // Not having indices is not necessarily an error, as the collection
  // could exist but just be empty
  if(m_sidreGroup->hasGroup(detail::COLLECTION_INDICES_NAME))
  {
    auto group = m_sidreGroup->getGroup(detail::COLLECTION_INDICES_NAME);
    indices.reserve(group->getNumViews());
    for(auto idx = group->getFirstValidViewIndex(); sidre::indexIsValid(idx);
        idx = group->getNextValidViewIndex(idx))
    {
      auto view = group->getView(idx);
      if(view->getTypeID() == axom::sidre::CHAR8_STR_ID)
      {
        std::string string_idx = view->getString();
        VariantKey key = string_idx;
        if(trimAbsolute)
        {
          // If the index is full/absolute, we only care about the last segment of it
          string_idx = removeBeforeDelimiter(string_idx);
          // The basename might be an integer, so check and convert accordingly
          int idx_as_int;
          if(checkedConvertToInt(string_idx, idx_as_int))
          {
            key = idx_as_int;
          }
          else
          {
            key = string_idx;
          }
        }
        indices.push_back(key);
      }
      else
      {
        indices.push_back(view->getData<int>());
      }
    }
  }
  return indices;
}

std::vector<std::pair<std::string, std::string>>
Container::collectionIndicesWithPaths(const std::string& name) const
{
  std::vector<std::pair<std::string, std::string>> result;
  for(const auto& indexLabel : collectionIndices(false))
  {
    auto stringLabel = detail::indexToString(indexLabel);
    // Since the index is absolute, we only care about the last segment of it
    // But since it's an absolute path then it gets used as the fullPath
    // which is used by the Reader to search in the input file
    const auto baseName = removeBeforeDelimiter(stringLabel);
    const auto fullPath = appendPrefix(stringLabel, name);
    // e.g. fullPath could be foo/1/bar for field "bar" at index 1 of array "foo"
    result.push_back({baseName, fullPath});
  }
  return result;
}

Verifiable<Container>& Container::addBoolArray(const std::string& name,
                                               const std::string& description)
{
  return addPrimitiveArray<bool>(name, description);
}

Verifiable<Container>& Container::addIntArray(const std::string& name,
                                              const std::string& description)
{
  return addPrimitiveArray<int>(name, description);
}

Verifiable<Container>& Container::addDoubleArray(const std::string& name,
                                                 const std::string& description)
{
  return addPrimitiveArray<double>(name, description);
}

Verifiable<Container>& Container::addStringArray(const std::string& name,
                                                 const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description);
}

template <typename Key>
Container& Container::addStructCollection(const std::string& name,
                                          const std::string& description)
{
  auto& container =
    addContainer(appendPrefix(name, detail::COLLECTION_GROUP_NAME), description);
  for(Container& sub_container : m_nested_aggregates)
  {
    container.m_nested_aggregates.push_back(
      sub_container.addStructCollection<Key>(name, description));
  }
  if(isStructCollection())
  {
    // Iterate over each element and forward the call to addPrimitiveArray
    for(const auto& indexPath : collectionIndicesWithPaths(name))
    {
      container.m_nested_aggregates.push_back(
        getContainer(indexPath.first).addStructCollection<Key>(name, description));
    }
    markAsStructCollection(*container.m_sidreGroup);
  }
  else
  {
    std::vector<Key> indices;
    std::string fullName = appendPrefix(m_name, name);
    fullName = removeAllInstances(fullName, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(fullName, m_unexpectedNames);
    const auto result = m_reader.getIndices(fullName, indices);
    if(result == ReaderResult::Success)
    {
      container.addIndicesGroup(indices, description);
    }
    markRetrievalStatus(*container.m_sidreGroup, result);
    markAsStructCollection(*container.m_sidreGroup);
  }
  return container;
}

Container& Container::addStructArray(const std::string& name,
                                     const std::string& description)
{
  return addStructCollection<int>(name, description);
}

Verifiable<Container>& Container::addBoolDictionary(const std::string& name,
                                                    const std::string& description)
{
  return addPrimitiveArray<bool>(name, description, true);
}

Verifiable<Container>& Container::addIntDictionary(const std::string& name,
                                                   const std::string& description)
{
  return addPrimitiveArray<int>(name, description, true);
}

Verifiable<Container>& Container::addDoubleDictionary(const std::string& name,
                                                      const std::string& description)
{
  return addPrimitiveArray<double>(name, description, true);
}

Verifiable<Container>& Container::addStringDictionary(const std::string& name,
                                                      const std::string& description)
{
  return addPrimitiveArray<std::string>(name, description, true);
}

Container& Container::addStructDictionary(const std::string& name,
                                          const std::string& description)
{
  return addStructCollection<VariantKey>(name, description);
}

axom::sidre::Group* Container::createSidreGroup(const std::string& name,
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

Field& Container::addField(axom::sidre::Group* sidreGroup,
                           axom::sidre::DataTypeId type,
                           const std::string& fullName,
                           const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currContainer = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Containers (if not present) before adding the field
    currContainer = &addContainer(name.substr(0, found));
  }
  const auto& emplace_result = currContainer->m_fieldChildren.emplace(
    fullName,
    cpp11_compat::make_unique<Field>(sidreGroup,
                                     m_sidreRootGroup,
                                     type,
                                     m_docEnabled));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

Function& Container::addFunctionInternal(axom::sidre::Group* sidreGroup,
                                         FunctionVariant&& func,
                                         const std::string& fullName,
                                         const std::string& name)
{
  const size_t found = name.find_last_of("/");
  auto currContainer = this;
  if(found != std::string::npos)
  {
    // This will add any intermediate Containers (if not present) before adding the field
    currContainer = &addContainer(name.substr(0, found));
  }
  const auto& emplace_result = currContainer->m_functionChildren.emplace(
    fullName,
    cpp11_compat::make_unique<Function>(sidreGroup,
                                        m_sidreRootGroup,
                                        std::move(func),
                                        m_docEnabled));
  // emplace_result is a pair whose first element is an iterator to the inserted element
  return *(emplace_result.first->second);
}

template <typename T, typename SFINAE>
VerifiableScalar& Container::addPrimitive(const std::string& name,
                                          const std::string& description,
                                          bool forArray,
                                          T val,
                                          const std::string& pathOverride)
{
  if(isStructCollection() || !m_nested_aggregates.empty())
  {
    // If it has indices, we're adding a primitive field to an array
    // of structs, so we need to iterate over the subcontainers
    // corresponding to elements of the array
    std::vector<std::reference_wrapper<VerifiableScalar>> fields;
    for(Container& container : m_nested_aggregates)
    {
      fields.push_back(
        container.addPrimitive<T>(name, description, forArray, val));
    }
    if(isStructCollection())
    {
      for(const auto& indexPath : collectionIndicesWithPaths(name))
      {
        // Add a primitive to an array element (which is a struct)
        fields.push_back(
          getContainer(indexPath.first)
            .addPrimitive<T>(name, description, forArray, val, indexPath.second));
      }
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
      fmt::format("Failed to create Sidre group with name '{0}'", fullName));
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    auto typeId = addPrimitiveHelper(sidreGroup, lookupPath, forArray, val);
    return addField(sidreGroup, typeId, fullName, name);
  }
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<bool>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  bool val)
{
  const auto result = m_reader.getBool(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val ? int8(1) : int8(0));
  }
  markRetrievalStatus(*sidreGroup, result);
  return axom::sidre::DataTypeId::INT8_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<int>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  int val)
{
  const auto result = m_reader.getInt(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val);
  }
  markRetrievalStatus(*sidreGroup, result);
  return axom::sidre::DataTypeId::INT_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<double>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  double val)
{
  const auto result = m_reader.getDouble(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewScalar("value", val);
  }
  markRetrievalStatus(*sidreGroup, result);
  return axom::sidre::DataTypeId::DOUBLE_ID;
}

template <>
axom::sidre::DataTypeId Container::addPrimitiveHelper<std::string>(
  axom::sidre::Group* sidreGroup,
  const std::string& lookupPath,
  bool forArray,
  std::string val)
{
  const auto result = m_reader.getString(lookupPath, val);
  if(forArray || result == ReaderResult::Success)
  {
    sidreGroup->createViewString("value", val);
  }
  markRetrievalStatus(*sidreGroup, result);
  return axom::sidre::DataTypeId::CHAR8_STR_ID;
}

namespace detail
{
/*!
  *****************************************************************************
  * \brief Adds the contents of an array to the container
  * 
  * \return The keys that were added
  *****************************************************************************
  */
template <typename T>
std::vector<VariantKey> registerCollection(Container& container,
                                           const std::unordered_map<int, T>& collection)
{
  std::vector<VariantKey> result;
  for(const auto& entry : collection)
  {
    result.push_back(entry.first);
    container.addPrimitive(std::to_string(entry.first), "", true, entry.second);
  }
  return result;
}

/*!
  *****************************************************************************
  * \brief Adds the contents of a dict to the container
  * 
  * \return The keys that were added
  *****************************************************************************
  */
template <typename T>
std::vector<VariantKey> registerCollection(
  Container& container,
  const std::unordered_map<VariantKey, T>& collection)
{
  std::vector<VariantKey> result;
  for(const auto& entry : collection)
  {
    result.push_back(entry.first);
    auto string_key = indexToString(entry.first);
    SLIC_ERROR_IF(
      string_key.find('/') != std::string::npos,
      fmt::format("[Inlet] Dictionary key '{0}' contains illegal character '/'",
                  string_key));
    SLIC_ERROR_IF(string_key.empty(),
                  "[Inlet] Dictionary key cannot be the empty string");
    container.addPrimitive(string_key, "", true, entry.second);
  }
  return result;
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
   * \brief Finalizes the creation of a collection
   * \param [inout] container The container to add the collection to
   * \param [in] reader The Reader object to read the collection from
   * \param [in] lookupPath The path within the input file to the collection
   * 
   * \return The keys from the collection that was added
   *****************************************************************************
   */
  static std::vector<VariantKey> add(Container& container,
                                     Reader& reader,
                                     const std::string& lookupPath)
  {
    std::unordered_map<Key, bool> map;
    // Failure to retrieve a map is not necessarily an error
    const auto result = reader.getBoolMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, int>
{
  static std::vector<VariantKey> add(Container& container,
                                     Reader& reader,
                                     const std::string& lookupPath)
  {
    std::unordered_map<Key, int> map;
    const auto result = reader.getIntMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, double>
{
  static std::vector<VariantKey> add(Container& container,
                                     Reader& reader,
                                     const std::string& lookupPath)
  {
    std::unordered_map<Key, double> map;
    const auto result = reader.getDoubleMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

template <typename Key>
struct PrimitiveArrayHelper<Key, std::string>
{
  static std::vector<VariantKey> add(Container& container,
                                     Reader& reader,
                                     const std::string& lookupPath)
  {
    std::unordered_map<Key, std::string> map;
    const auto result = reader.getStringMap(lookupPath, map);
    markRetrievalStatus(*container.sidreGroup(), result);
    return registerCollection(container, map);
  }
};

void addIndexViewToGroup(sidre::Group& group, const int& index)
{
  group.createViewScalar("", index);
}

void addIndexViewToGroup(sidre::Group& group, const std::string& index)
{
  group.createViewString("", index);
}

void addIndexViewToGroup(sidre::Group& group, const VariantKey& index)
{
  if(index.type() == InletType::String)
  {
    addIndexViewToGroup(group, static_cast<std::string>(index));
  }
  else
  {
    addIndexViewToGroup(group, static_cast<int>(index));
  }
}

void updateUnexpectedNames(const std::string& accessedName,
                           std::unordered_set<std::string>& unexpectedNames)
{
  for(auto iter = unexpectedNames.begin(); iter != unexpectedNames.end();)
  {
    // Check if the possibly unexpected name is a substring of the accessed name,
    // if it is, then it gets marked as expected via removal
    if(accessedName.find(*iter) != std::string::npos)
    {
      iter = unexpectedNames.erase(iter);
    }
    else
    {
      ++iter;
    }
  }
}

}  // end namespace detail

template <typename Key>
void Container::addIndicesGroup(const std::vector<Key>& indices,
                                const std::string& description)
{
  sidre::Group* indices_group =
    m_sidreGroup->createGroup(detail::COLLECTION_INDICES_NAME,
                              /* list_format = */ true);
  // For each index, add a container whose name is its index
  // Schema for struct is defined using the returned container
  for(const auto& idx : indices)
  {
    const std::string string_idx =
      removeBeforeDelimiter(detail::indexToString(idx));
    addContainer(string_idx, description);
    std::string absolute = appendPrefix(m_name, detail::indexToString(idx));
    absolute = removeAllInstances(absolute, detail::COLLECTION_GROUP_NAME + "/");
    detail::addIndexViewToGroup(*indices_group, absolute);
  }
}

template <typename T, typename SFINAE>
Verifiable<Container>& Container::addPrimitiveArray(const std::string& name,
                                                    const std::string& description,
                                                    const bool isDict,
                                                    const std::string& pathOverride)
{
  if(isStructCollection() || !m_nested_aggregates.empty())
  {
    // Adding an array of primitive field to an array of structs
    std::vector<std::reference_wrapper<Verifiable>> containers;
    for(Container& container : m_nested_aggregates)
    {
      containers.push_back(
        container.addPrimitiveArray<T>(name, description, isDict));
    }
    if(isStructCollection())
    {
      // Iterate over each element and forward the call to addPrimitiveArray
      for(const auto& indexPath : collectionIndicesWithPaths(name))
      {
        containers.push_back(
          getContainer(indexPath.first)
            .addPrimitiveArray<T>(name, description, isDict, indexPath.second));
      }
    }

    m_aggregate_containers.emplace_back(std::move(containers));

    // Remove when C++17 is available
    return m_aggregate_containers.back();
  }
  else
  {
    // "base case", create a container for the field and fill it in with the helper
    auto& container =
      addContainer(appendPrefix(name, detail::COLLECTION_GROUP_NAME),
                   description);
    const std::string& fullName = appendPrefix(m_name, name);
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    std::vector<VariantKey> indices;
    if(isDict)
    {
      indices = detail::PrimitiveArrayHelper<VariantKey, T>::add(container,
                                                                 m_reader,
                                                                 lookupPath);
    }
    else
    {
      indices =
        detail::PrimitiveArrayHelper<int, T>::add(container, m_reader, lookupPath);
    }
    // Copy the indices to the datastore to keep track of integer vs. string indices
    if(!indices.empty())
    {
      container.addIndicesGroup(indices, description);
    }
    return container;
  }
}

Verifiable<Function>& Container::addFunction(const std::string& name,
                                             const FunctionTag ret_type,
                                             const std::vector<FunctionTag>& arg_types,
                                             const std::string& description,
                                             const std::string& pathOverride)
{
  if(isStructCollection() || !m_nested_aggregates.empty())
  {
    // If it has indices, we're adding a primitive field to an array
    // of structs, so we need to iterate over the subcontainers
    // corresponding to elements of the array
    std::vector<std::reference_wrapper<Verifiable<Function>>> funcs;
    for(Container& container : m_nested_aggregates)
    {
      funcs.push_back(
        container.addFunction(name, ret_type, arg_types, description));
    }
    if(isStructCollection())
    {
      for(const auto& indexPath : collectionIndicesWithPaths(name))
      {
        // Add a primitive to an array element (which is a struct)
        funcs.push_back(
          getContainer(indexPath.first)
            .addFunction(name, ret_type, arg_types, description, indexPath.second));
      }
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
      fmt::format("Failed to create Sidre group with name '{0}'", fullName));
    // If a pathOverride is specified, needed when Inlet-internal groups
    // are part of fullName
    std::string lookupPath = (pathOverride.empty()) ? fullName : pathOverride;
    lookupPath =
      removeAllInstances(lookupPath, detail::COLLECTION_GROUP_NAME + "/");
    detail::updateUnexpectedNames(lookupPath, m_unexpectedNames);
    auto func = m_reader.getFunction(lookupPath, ret_type, arg_types);
    return addFunctionInternal(sidreGroup, std::move(func), fullName, name);
  }
}

Proxy Container::operator[](const std::string& name) const
{
  const bool has_container = hasContainer(name);
  const bool has_field = hasField(name);
  const bool has_func = hasFunction(name);

  // Ambiguous case - both a container and field exist with the same name
  if((has_container && has_field) || (has_field && has_func) ||
     (has_container && has_func))
  {
    const std::string msg = fmt::format(
      "[Inlet] Ambiguous lookup - more than one of a container/field/function "
      "with name '{0}' exist",
      name);
    SLIC_ERROR(msg);
    return Proxy();
  }

  else if(has_container)
  {
    return Proxy(getContainer(name));
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
    std::string msg = fmt::format(
      "[Inlet] No container, field, or function with name '{0}' exists",
      name);
    SLIC_ERROR(msg);
    return Proxy();
  }
}

Container& Container::required(bool isRequired)
{
  // If it's a struct collection we set the individual fields as required,
  // and also the collection container itself, as the user would expect that marking
  // a struct collection as required means that it is non-empty
  if(isStructCollection())
  {
    forEachCollectionElement(
      [isRequired](Container& container) { container.required(isRequired); });
  }

  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Container specific Sidre Datastore Group not set");
  setRequired(*m_sidreGroup, *m_sidreRootGroup, isRequired);
  return *this;
}

bool Container::isRequired() const
{
  if(isStructCollection())
  {
    bool result = false;
    forEachCollectionElement([&result](Container& container) {
      if(container.isRequired())
      {
        result = true;
      }
    });
    return result;
  }

  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Container specific Sidre Datastore Group not set");
  return checkIfRequired(*m_sidreGroup, *m_sidreRootGroup);
}

Container& Container::registerVerifier(std::function<bool(const Container&)> lambda)
{
  if(isStructCollection())
  {
    forEachCollectionElement(
      [&lambda](Container& container) { container.registerVerifier(lambda); });
  }
  else
  {
    SLIC_WARNING_IF(m_verifier,
                    fmt::format("[Inlet] Verifier for Container "
                                "already set: {0}",
                                m_name));
    m_verifier = lambda;
  }
  return *this;
}

bool Container::verify() const
{
  // Whether the calling container has any "truthy" subcontainers, fields, or functions
  // If the name is empty then we're the global (root) container, which we always
  // consider to be defined
  const bool this_container_defined = static_cast<bool>(*this) || m_name.empty();

  // If this container was required, make sure something was defined in it
  bool verified =
    verifyRequired(*m_sidreGroup, this_container_defined, "Container");

  // Verify this Container if a lambda was configured
  if(this_container_defined && m_verifier && !m_verifier(*this))
  {
    verified = false;
    SLIC_WARNING(
      fmt::format("[Inlet] Container failed verification: {0}", m_name));
  }

  // Checking the child objects is not needed if the container is empty
  if(this_container_defined)
  {
    // Verify the child Fields of this Container
    for(const auto& field : m_fieldChildren)
    {
      verified = verified && field.second->verify();
    }
    // Verify the child Containers of this Container
    for(const auto& container : m_containerChildren)
    {
      verified = verified && container.second->verify();
    }

    // Verify the child Functions of this Container
    for(const auto& function : m_functionChildren)
    {
      verified = verified && function.second->verify();
    }
  }
  // If this has a collection group, it always needs to be verified, as annotations
  // may have been applied to the collection group and not the calling group
  else if(hasContainer(detail::COLLECTION_GROUP_NAME))
  {
    verified = verified && getContainer(detail::COLLECTION_GROUP_NAME).verify();
  }

  return verified;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Container>> Container::*
Container::getChildren()
{
  return &Container::m_containerChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Field>> Container::*
Container::getChildren()
{
  return &Container::m_fieldChildren;
}

template <>
std::unordered_map<std::string, std::unique_ptr<Function>> Container::*
Container::getChildren()
{
  return &Container::m_functionChildren;
}

template <typename T>
bool Container::hasChild(const std::string& childName) const
{
  const auto& children = this->*getChildren<T>();
  return children.find(appendPrefix(m_name, childName)) != children.end();
}

template <typename T>
T* Container::getChildInternal(const std::string& childName) const
{
  std::string name = childName;
  size_t found = name.find("/");
  auto currContainer = this;

  while(found != std::string::npos)
  {
    const std::string& currName = name.substr(0, found);
    if(currContainer->hasChild<Container>(currName))
    {
      currContainer = currContainer->m_containerChildren
                        .at(appendPrefix(currContainer->m_name, currName))
                        .get();
    }
    else
    {
      return nullptr;
    }
    name = name.substr(found + 1);
    found = name.find("/");
  }

  if(currContainer->hasChild<T>(name))
  {
    const auto& children = currContainer->*getChildren<T>();
    return children.at(appendPrefix(currContainer->m_name, name)).get();
  }
  return nullptr;
}

bool Container::hasContainer(const std::string& containerName) const
{
  return static_cast<bool>(getChildInternal<Container>(containerName));
}

bool Container::hasField(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Field>(fieldName));
}

bool Container::hasFunction(const std::string& fieldName) const
{
  return static_cast<bool>(getChildInternal<Function>(fieldName));
}

Container& Container::getContainer(const std::string& containerName) const
{
  auto container = getChildInternal<Container>(containerName);
  if(!container)
  {
    SLIC_ERROR(fmt::format("[Inlet] Container not found: {0}", containerName));
  }
  return *container;
}

Field& Container::getField(const std::string& fieldName) const
{
  auto field = getChildInternal<Field>(fieldName);
  if(!field)
  {
    SLIC_ERROR(fmt::format("[Inlet] Field not found: {0}", fieldName));
  }
  return *field;
}

Function& Container::getFunction(const std::string& funcName) const
{
  auto func = getChildInternal<Function>(funcName);
  if(!func)
  {
    SLIC_ERROR(fmt::format("[Inlet] Function not found: {0}", funcName));
  }
  return *func;
}

std::string Container::name() const { return m_name; }

bool Container::contains(const std::string& name) const
{
  if(auto container = getChildInternal<Container>(name))
  {
    // call operator bool on the container itself
    return static_cast<bool>(*container);
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

Container::operator bool() const
{
  // Check if any of its child containers are nontrivial
  const bool has_containers =
    std::any_of(m_containerChildren.begin(),
                m_containerChildren.end(),
                [](const decltype(m_containerChildren)::value_type& entry) {
                  return static_cast<bool>(*entry.second);
                });

  const bool has_fields =
    std::any_of(m_fieldChildren.begin(),
                m_fieldChildren.end(),
                [](const decltype(m_fieldChildren)::value_type& entry) {
                  return static_cast<bool>(*entry.second);
                });

  const bool has_functions =
    std::any_of(m_functionChildren.begin(),
                m_functionChildren.end(),
                [](const decltype(m_functionChildren)::value_type& entry) {
                  return static_cast<bool>(*entry.second);
                });

  return has_containers || has_fields || has_functions;
}

const std::unordered_map<std::string, std::unique_ptr<Container>>&
Container::getChildContainers() const
{
  return m_containerChildren;
}

const std::unordered_map<std::string, std::unique_ptr<Field>>&
Container::getChildFields() const
{
  return m_fieldChildren;
}

}  // end namespace inlet
}  // end namespace axom
