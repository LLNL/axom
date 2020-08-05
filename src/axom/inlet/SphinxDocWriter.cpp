// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxDocWriter.cpp
 *
 * \brief This file contains the class implementation of the SphinxDocWriter.
 *******************************************************************************
 */

#include "axom/inlet/SphinxDocWriter.hpp"
#include <iostream>
#include <assert.h>
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
 
SphinxDocWriter::SphinxDocWriter(const std::string& fileName, axom::sidre::Group* root) {
  SLIC_WARNING_IF(root == nullptr, " Inlet's Sidre root group was not defined"
                                   " while writing documentation");
  m_sidreRootGroup = root;
  m_fileName = fileName;
}

void SphinxDocWriter::writeDocumentation() {
  if (m_sidreRootGroup->getName() == "") {
    writeTitle("Input Deck Options");
  }  else {
    writeTitle(m_sidreRootGroup->getName());
  }
  m_oss << ".. |uncheck|    unicode:: U+2610 .. UNCHECKED BOX" << std::endl;
  m_oss << ".. |check|      unicode:: U+2611 .. CHECKED BOX" << std::endl;
  writeDocumentationHelper(m_sidreRootGroup);
  writeAllTables();
  m_outFile.open(m_fileName);
  m_outFile << m_oss.str();
  m_outFile.close();
}

void SphinxDocWriter::writeDocumentationHelper(axom::sidre::Group* sidreGroup) {
  SLIC_WARNING_IF(!sidreGroup, "Root was nullptr");
  axom::sidre::IndexType i = sidreGroup->getFirstValidGroupIndex();

  // Special case for valid string values since it's a Group whereas other
  // valid values and Field properties are all stored as Views
  if (axom::sidre::indexIsValid(i) && sidreGroup->getGroupName(i) == "validStringValues") {
    i = axom::sidre::InvalidIndex;
  } 

  // Case 1: the current group is a Field
  if (sidreGroup != m_sidreRootGroup && i == axom::sidre::InvalidIndex) { 
    extractFieldMetadata(sidreGroup);
  } 

  // Case 2: Current root corresponds to an Inlet::Table
  if (i != axom::sidre::InvalidIndex) {
    m_inletTablePathNames.push_back(sidreGroup->getPathName());
    m_rstTables[sidreGroup->getPathName()] = TableData();
    m_rstTables[sidreGroup->getPathName()].tableName = sidreGroup->getName();
    if (sidreGroup->getName() != "" && sidreGroup->hasView("description")) {
      m_rstTables[sidreGroup->getPathName()].description = sidreGroup->getView("description")->getString();
    } 
  }

  while (i != axom::sidre::InvalidIndex) {
    writeDocumentationHelper(sidreGroup->getGroup(i));
    i = sidreGroup->getNextValidGroupIndex(i);
  }
}

void SphinxDocWriter::writeTitle(const std::string& title) {
  if (title != "") {
    std::string equals = std::string(title.length(), '=');
    m_oss << equals << "\n" << title << "\n" << equals << "\n";
  }
}

void SphinxDocWriter::writeSubtitle(const std::string& sub) {
  if (sub != "") {
    std::string dashes = std::string(sub.length(), '-');
    m_oss << "\n" << dashes << "\n" << sub << "\n" << dashes << "\n\n";
  }
}

void SphinxDocWriter::writeTable(const std::string& title, 
                                 const std::vector<std::vector<std::string>>& rstTable) {
  SLIC_WARNING_IF(rstTable.size() <= 1, "Vector for corresponding rst table must be nonempty");
  std::string result = ".. list-table:: " + title;
  result += "\n   :widths: 25 25 25 25 25\n   :header-rows: 1\n   :stub-columns: 1\n\n";
  for (unsigned int i = 0; i < rstTable.size(); ++i) {
    result += "   * - ";
    for (unsigned int j = 0; j < rstTable[i].size(); ++j) {
      if (j != 0) {
        result += "     - ";
      }
      result += rstTable[i][j] + "\n";
    }
  }
  m_oss << result;
}

void SphinxDocWriter::writeAllTables() {
  for (std::string& pathName : m_inletTablePathNames) {
    writeSubtitle(m_rstTables[pathName].tableName);
    if (m_rstTables[pathName].description != "") {
      m_oss << "Description: " << m_rstTables[pathName].description 
            << std::endl << std::endl;
    }
    if (m_rstTables[pathName].rstTable.size() > 1) {
      writeTable("Fields", m_rstTables[pathName].rstTable);
    }
  }
}

std::string SphinxDocWriter::getDefaultValueAsString(axom::sidre::View* view) {
  axom::sidre::TypeID type = view->getTypeID(); 
  if (type == axom::sidre::TypeID::INT8_ID) {
    int8 val = view->getData();
    return val ? "True" : "False";
  } else if (type == axom::sidre::TypeID::INT_ID) {
    int val = view->getData();
    return std::to_string(val);
  } else if (type == axom::sidre::TypeID::DOUBLE_ID) {
    double val = view->getData();
    return std::to_string(val);
  } 
  return view->getString();
}

std::string SphinxDocWriter::getRangeAsString(axom::sidre::View* view) {
  std::ostringstream oss;
  oss.precision(3);
  oss << std::scientific;

  axom::sidre::TypeID type = view->getTypeID();
  if (type == axom::sidre::INT_ID) {
    int* range = view->getArray();
    oss << range[0] << " to " << range[1];
  } else {
    double* range = view->getArray();
    oss << range[0] << " to " << range[1];
  }
  return oss.str();
}

std::string SphinxDocWriter::getValidValuesAsString(axom::sidre::View* view) {
  int* range = view->getArray();
  size_t size = view->getBuffer()->getNumElements();
  std::string result = "";
  for (size_t i = 0; i < size; ++i) {
    if (i == size-1) {
      result += std::to_string(range[i]);      
    } else {
      result += std::to_string(range[i]) + ", ";
    }
  }
  return result;
}

std::string SphinxDocWriter::getValidStringValues(axom::sidre::Group* sidreGroup) {
  auto idx = sidreGroup->getFirstValidViewIndex();
  std::string validValues = "";
  while (axom::sidre::indexIsValid(idx)) {
    validValues += std::string(sidreGroup->getView(idx)->getString());
    idx = sidreGroup->getNextValidViewIndex(idx);
    if (axom::sidre::indexIsValid(idx)) {
      validValues += ", ";
    }
  }
  return validValues;
}

void SphinxDocWriter::extractFieldMetadata(axom::sidre::Group* sidreGroup) {
  TableData& currentTable = m_rstTables[sidreGroup->getParent()->getPathName()];
  std::vector<std::string> fieldAttributes(5, "");
  fieldAttributes.resize(5);

  fieldAttributes[0] = sidreGroup->getName();

  if (sidreGroup->hasView("description")) {
    fieldAttributes[1] = std::string(sidreGroup->getView("description")->getString());
  } 

  if (sidreGroup->hasView("defaultValue")) { 
    fieldAttributes[2] = getDefaultValueAsString(sidreGroup->getView("defaultValue"));
  }

  if (sidreGroup->hasView("range")) {
    fieldAttributes[3] = getRangeAsString(sidreGroup->getView("range"));
  } else if (sidreGroup->hasView("validValues")) {
    fieldAttributes[3] = getValidValuesAsString(sidreGroup->getView("validValues"));
  } else if (sidreGroup->hasGroup("validStringValues")) {
    fieldAttributes[3] = getValidStringValues(sidreGroup->getGroup("validStringValues"));
  }
  
  if (sidreGroup->hasView("required")) {
    int8 required = sidreGroup->getView("required")->getData();
    fieldAttributes[4] = required ? "|check|" : "|uncheck|";
  } else {
    fieldAttributes[4] = "|uncheck|";
  }

  currentTable.rstTable.push_back(fieldAttributes);
}

}
}
