// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxWriter.cpp
 *
 * \brief This file contains the class implementation of the SphinxWriter.
 *******************************************************************************
 */

#include "axom/inlet/SphinxWriter.hpp"

#include <iostream>

#include "axom/slic.hpp"
#include "axom/inlet/Container.hpp"

namespace axom
{
namespace inlet
{
namespace detail
{
/**
 * \brief Determines whether a container is trivial (contains no fields/functions in its subtree)
 * 
 * \param [in] container The container to evaluate
 */
bool isTrivial(const Container& container)
{
  if(!container.getChildFields().empty() ||
     !container.getChildFunctions().empty())
  {
    return false;
  }
  using value_type =
    std::decay<decltype(container.getChildContainers())>::type::value_type;
  return std::all_of(
    container.getChildContainers().begin(),
    container.getChildContainers().end(),
    [](const value_type& entry) { return isTrivial(*entry.second); });
}

/**
 * \brief Converts an enumeration to its underlying type
 * \param [in] e The enumeration value to convert
 * This function should be removed once C++23 is available
 * \see https://en.cppreference.com/w/cpp/utility/to_underlying
 */
template <typename E>
constexpr typename std::underlying_type<E>::type to_underlying(const E e)
{
  return static_cast<typename std::underlying_type<E>::type>(e);
}

}  // namespace detail

SphinxWriter::SphinxWriter(const std::string& fileName,
                           const std::string& title,
                           const Style style)
  : m_fieldColLabels({"Field Name",
                      "Description",
                      "Default Value",
                      "Range/Valid Values",
                      "Required"})
  , m_functionColLabels(
      {"Function Name", "Description", "Signature", "Required"})
  , m_style(style)
{
  m_fileName = fileName;
  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX\n";
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX\n\n";
  if(title.empty())
  {
    writeTitle("Input file Options");
  }
  else
  {
    writeTitle(title);
  }
}

void SphinxWriter::documentContainer(const Container& container)
{
  const auto sidreGroup = container.sidreGroup();
  const std::string pathName = sidreGroup->getPathName();
  std::string containerName = sidreGroup->getName();
  bool isSelectedElement = false;

  // If the container is empty, ignore it
  if(detail::isTrivial(container))
  {
    return;
  }

  // Replace the "implementation-defined" name with something a bit more readable
  if(isCollectionGroup(containerName))
  {
    containerName = "Collection contents:";
  }

  // If we've gotten to this point and are an element of an array/dict,
  // mark it as the selected element
  if(sidreGroup->getParent()->getName() == detail::COLLECTION_GROUP_NAME)
  {
    // The collection that this Container is a part of
    const std::string collectionName =
      sidreGroup->getParent()->getParent()->getPathName();
    isSelectedElement = true;
  }

  m_inletContainerPathNames.push_back(pathName);
  auto& currContainer =
    m_rstTables
      .emplace(pathName, ContainerData {m_fieldColLabels, m_functionColLabels})
      .first->second;
  currContainer.containerName = containerName;
  currContainer.isSelectedElement = isSelectedElement;
  if(containerName != "" && sidreGroup->hasView("description"))
  {
    currContainer.description = sidreGroup->getView("description")->getString();
  }

  // Bail out if this is a collection of primitives -
  // it doesn't make sense to document each individual primitive
  // FIXME: Implement something analogous to the logic for struct collections that displays the
  // "schema" for the primitive elements of the collection
  if(isCollectionGroup(container.name()) &&
     !sidreGroup->hasView(detail::STRUCT_COLLECTION_FLAG))
  {
    return;
  }

  for(const auto& field_entry : container.getChildFields())
  {
    extractFieldMetadata(field_entry.second->sidreGroup(), currContainer);
  }

  for(const auto& function_entry : container.getChildFunctions())
  {
    extractFunctionMetadata(function_entry.second->sidreGroup(), currContainer);
  }

  // If it's not a collection,  we need to record the child containers
  if(!isCollectionGroup(container.name()))
  {
    for(const auto& container_entry : container.getChildContainers())
    {
      if(!isCollectionGroup(container_entry.first))
      {
        const auto name = removeBeforeDelimiter(container_entry.first);
        std::string description;
        const auto group = container_entry.second->sidreGroup();
        const static auto collectionDescriptionPath =
          appendPrefix(detail::COLLECTION_GROUP_NAME, "description");
        if(group->hasView("description"))
        {
          description = group->getView("description")->getString();
        }
        else if(group->hasView(collectionDescriptionPath))
        {
          // If the label is applied to the collection group itself, we can use that instead
          description = group->getView(collectionDescriptionPath)->getString();
        }
        currContainer.childContainers.push_back(
          {std::move(name),
           std::move(description),
           detail::isTrivial(*container_entry.second)});
      }
    }
  }
}

void SphinxWriter::finalize()
{
  if(m_style == Style::Nested)
  {
    writeNestedTables();
  }
  else
  {
    writeAllTables();
  }
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxWriter::writeTitle(const std::string& title)
{
  if(title != "")
  {
    std::string equals = std::string(title.length(), '=');
    m_oss << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxWriter::writeSubtitle(const std::string& sub)
{
  if(sub != "")
  {
    std::string dashes = std::string(sub.length(), '-');
    m_oss << "\n" << dashes << "\n" << sub << "\n" << dashes << "\n\n";
  }
}

void SphinxWriter::writeTable(const std::string& title,
                              const std::vector<std::vector<std::string>>& rstTable)
{
  SLIC_WARNING_IF(
    rstTable.size() <= 1,
    "[Inlet] Vector for corresponding rst table must be nonempty");
  std::string result = ".. list-table:: " + title;
  std::string widths = ":widths:";
  // This would be easier with an iterator adaptor like back_inserter but for
  // concatenation
  for(std::size_t i = 0u; i < rstTable.front().size(); i++)
  {
    widths += " 25";
  }
  result += "\n   " + widths + "\n";
  result += "   :header-rows: 1\n   :stub-columns: 1\n\n";
  for(unsigned int i = 0; i < rstTable.size(); ++i)
  {
    result += "   * - ";
    for(unsigned int j = 0; j < rstTable[i].size(); ++j)
    {
      if(j != 0)
      {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  m_oss << result;
}

void SphinxWriter::writeAllTables()
{
  for(std::string& pathName : m_inletContainerPathNames)
  {
    auto& currContainer = m_rstTables.at(pathName);
    // If we're displaying a selected element, the title and description
    // will already have been printed
    if(currContainer.isSelectedElement)
    {
      m_oss << "The input schema defines a collection of this container.\n";
      m_oss << "For brevity, only one instance is displayed here.\n\n";
    }
    else
    {
      writeSubtitle(currContainer.containerName);
      if(currContainer.description != "")
      {
        m_oss << "Description: " << currContainer.description << "\n\n";
      }
    }
    if(currContainer.fieldTable.size() > 1)
    {
      writeTable("Fields", currContainer.fieldTable);
    }
    if(currContainer.functionTable.size() > 1)
    {
      writeTable("Functions", currContainer.functionTable);
    }
  }
}

void SphinxWriter::writeNestedTables()
{
  // Used to avoid duplicate hyperlinks
  int containerCount = 0;
  for(const auto& pathName : m_inletContainerPathNames)
  {
    const auto& currContainer = m_rstTables.at(pathName);

    if(!currContainer.isSelectedElement)
    {
      writeSubtitle(currContainer.containerName);
    }

    if(currContainer.description != "")
    {
      m_oss << currContainer.description << "\n\n";
    }

    const auto& fields = currContainer.fieldTable;
    const auto& functions = currContainer.functionTable;
    const auto& containers =
      currContainer.childContainers;  // Contains only names and descriptions

    if(fields.size() > 1 || functions.size() > 1 || !containers.empty())
    {
      m_oss << ".. list-table::\n";
      m_oss << "   :widths: 25 50\n";
      m_oss << "   :header-rows: 1\n";
      m_oss << "   :stub-columns: 1\n\n";
      m_oss << "   * - Name\n";
      m_oss << "     - Description\n";
    }

    // Writes a name + description to the "table of contents"
    // for each Container
    auto writeTOCEntry = [this,
                          containerCount](const std::vector<std::string>& data) {
      // FIXME: Overlapping link names? Need to use the full path or a hash??
      m_oss << fmt::format("   * - :ref:`{0}<{0}{1}>`\n", data[0], containerCount);
      m_oss << fmt::format("     - {0}\n", data[1]);
    };

    // FIXME: Would it be useful to alphabetize these??
    for(const auto& container : containers)
    {
      // We can set up a hyperlink to non-trivial tables because they will
      // be subsequently documented
      if(container.isTrivial)
      {
        m_oss << fmt::format("   * - {0}\n", container.name);
      }
      else
      {
        m_oss << fmt::format("   * - `{0}`_\n", container.name);
      }
      m_oss << fmt::format("     - {0}\n", container.description);
    }

    // Need to skip first element (header row), hence no range-based for loop
    std::for_each(fields.begin() + 1, fields.end(), writeTOCEntry);
    std::for_each(functions.begin() + 1, functions.end(), writeTOCEntry);

    m_oss << "\n\n";

    // Now, dump the full descriptions
    std::for_each(
      fields.begin() + 1,
      fields.end(),
      [this, containerCount](const std::vector<std::string>& fieldData) {
        // Insert a hyperlink target for the name
        m_oss << fmt::format(".. _{0}{1}:\n\n", fieldData[0], containerCount);
        m_oss << fmt::format("**{0}**\n\n", fieldData[0]);  // name
        // description - ideally this would be an extended description
        m_oss << fieldData[1] << "\n\n";

        // The default value
        if(!fieldData[2].empty())
        {
          m_oss << fmt::format("  - Default value: {0}\n", fieldData[2]);
        }

        // The valid values
        if(!fieldData[3].empty())
        {
          m_oss << fmt::format("  - Valid values: {0}\n", fieldData[3]);
        }

        // Whether it's required
        if(!fieldData[4].empty())
        {
          const bool isRequired = fieldData[3] == "|check|";
          m_oss << fmt::format("  - {0}\n", isRequired ? "Required" : "Optional");
        }
        m_oss << "\n\n";
      });

    std::for_each(
      functions.begin() + 1,
      functions.end(),
      [this, containerCount](const std::vector<std::string>& functionData) {
        // Insert a hyperlink target for the name
        m_oss << fmt::format(".. _{0}{1}:\n\n", functionData[0], containerCount);
        m_oss << fmt::format("**{0}**\n\n", functionData[0]);  // name
        // description - ideally this would be an extended description
        m_oss << functionData[1] << "\n\n";

        // The function's signature
        if(!functionData[2].empty())
        {
          m_oss << fmt::format("  - Signature: {0}\n", functionData[2]);
        }

        // Whether it's required
        if(!functionData[3].empty())
        {
          const bool isRequired = functionData[3] == "|check|";
          m_oss << fmt::format("  - {0}\n", isRequired ? "Required" : "Optional");
        }
        m_oss << "\n\n";
      });
    containerCount++;
  }
}

std::string SphinxWriter::getValueAsString(const axom::sidre::View* view)
{
  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::TypeID::INT8_ID)
  {
    int8 val = view->getData();
    return val ? "True" : "False";
  }
  else if(type == axom::sidre::TypeID::INT_ID)
  {
    int val = view->getData();
    return std::to_string(val);
  }
  else if(type == axom::sidre::TypeID::DOUBLE_ID)
  {
    double val = view->getData();
    return std::to_string(val);
  }
  return view->getString();
}

std::string SphinxWriter::getRangeAsString(const axom::sidre::View* view)
{
  std::ostringstream oss;
  oss.precision(3);
  oss << std::scientific;

  axom::sidre::TypeID type = view->getTypeID();
  if(type == axom::sidre::INT_ID)
  {
    const int* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  else
  {
    const double* range = view->getData();
    oss << range[0] << " to " << range[1];
  }
  return oss.str();
}

std::string SphinxWriter::getValidValuesAsString(const axom::sidre::View* view)
{
  const int* range = view->getData();
  size_t size = view->getBuffer()->getNumElements();
  std::string result = "";
  for(size_t i = 0; i < size; ++i)
  {
    if(i == size - 1)
    {
      result += std::to_string(range[i]);
    }
    else
    {
      result += std::to_string(range[i]) + ", ";
    }
  }
  return result;
}

std::string SphinxWriter::getValidStringValues(const axom::sidre::Group* sidreGroup)
{
  auto idx = sidreGroup->getFirstValidViewIndex();
  std::string validValues = "";
  while(axom::sidre::indexIsValid(idx))
  {
    validValues += std::string(sidreGroup->getView(idx)->getString());
    idx = sidreGroup->getNextValidViewIndex(idx);
    if(axom::sidre::indexIsValid(idx))
    {
      validValues += ", ";
    }
  }
  return validValues;
}

void SphinxWriter::extractFieldMetadata(const axom::sidre::Group* sidreGroup,
                                        ContainerData& currentContainer)
{
  std::vector<std::string> fieldAttributes(m_fieldColLabels.size());

  fieldAttributes[0] = sidreGroup->getName();

  if(sidreGroup->hasView("description"))
  {
    fieldAttributes[1] =
      std::string(sidreGroup->getView("description")->getString());
  }

  if(sidreGroup->hasView("defaultValue"))
  {
    fieldAttributes[2] = getValueAsString(sidreGroup->getView("defaultValue"));
  }

  if(sidreGroup->hasView("range"))
  {
    fieldAttributes[3] = getRangeAsString(sidreGroup->getView("range"));
  }
  else if(sidreGroup->hasView("validValues"))
  {
    fieldAttributes[3] =
      getValidValuesAsString(sidreGroup->getView("validValues"));
  }
  else if(sidreGroup->hasGroup("validStringValues"))
  {
    fieldAttributes[3] =
      getValidStringValues(sidreGroup->getGroup("validStringValues"));
  }

  if(sidreGroup->hasView("required"))
  {
    int8 required = sidreGroup->getView("required")->getData();
    fieldAttributes[4] = required ? "|check|" : "|uncheck|";
  }
  else
  {
    fieldAttributes[4] = "|uncheck|";
  }

  currentContainer.fieldTable.push_back(fieldAttributes);
}

std::string SphinxWriter::getSignatureAsString(const axom::sidre::Group* sidreGroup)
{
  using underlying = std::underlying_type<FunctionTag>::type;
  static const auto type_names = []() {
    std::unordered_map<underlying, std::string> result;
    result[detail::to_underlying(FunctionTag::Vector)] = "Vector";
    result[detail::to_underlying(FunctionTag::Double)] = "Double";
    result[detail::to_underlying(FunctionTag::Void)] = "Void";
    result[detail::to_underlying(FunctionTag::String)] = "String";
    return result;
  }();

  // View::getData<T> does not have a const version...
  const auto ret_type =
    static_cast<underlying>(sidreGroup->getView("return_type")->getData());

  const auto args_view = sidreGroup->getView("function_arguments");
  const underlying* arg_tags = args_view->getData();
  const int num_args = args_view->getNumElements();
  std::vector<std::string> arg_types(num_args);
  for(int i = 0; i < num_args; i++)
  {
    arg_types[i] = type_names.at(arg_tags[i]);
  }
  return fmt::format("{0}({1})",
                     type_names.at(ret_type),
                     fmt::join(arg_types, ", "));
}

void SphinxWriter::extractFunctionMetadata(const axom::sidre::Group* sidreGroup,
                                           ContainerData& currentContainer)
{
  std::vector<std::string> functionAttributes(m_functionColLabels.size());

  functionAttributes[0] = sidreGroup->getName();

  if(sidreGroup->hasView("description"))
  {
    functionAttributes[1] =
      std::string(sidreGroup->getView("description")->getString());
  }

  functionAttributes[2] = getSignatureAsString(sidreGroup);

  if(sidreGroup->hasView("required"))
  {
    int8 required = sidreGroup->getView("required")->getData();
    functionAttributes[3] = required ? "|check|" : "|uncheck|";
  }
  else
  {
    functionAttributes[3] = "|uncheck|";
  }

  currentContainer.functionTable.push_back(functionAttributes);
}

}  // namespace inlet
}  // namespace axom
