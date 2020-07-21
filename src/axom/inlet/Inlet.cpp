// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.cpp
 *
 * \brief This file contains the class implementation of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#include "axom/inlet/Inlet.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

std::shared_ptr<Table> Inlet::addTable(const std::string& name,
                                       const std::string& description)
{
  return std::make_shared<Table>(name, description, m_reader, m_sidreRootGroup, m_docEnabled);
}

std::shared_ptr<Field> Inlet::addBool(const std::string& name,
                                      const std::string& description)
{
  return m_globalTable->addBool(name, description);
}

std::shared_ptr<Field> Inlet::addDouble(const std::string& name,
                                        const std::string& description)
{
  return m_globalTable->addDouble(name, description);
}

std::shared_ptr<Field> Inlet::addInt(const std::string& name,
                                     const std::string& description)
{
  return m_globalTable->addInt(name, description);
}

std::shared_ptr<Field> Inlet::addString(const std::string& name,
                                        const std::string& description)
{
  return m_globalTable->addString(name, description);
}


//-------------------------------------------------
//   Get values out of the datastore
//-------------------------------------------------

axom::sidre::View* Inlet::baseGet(const std::string& name)
{
  SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore Group not set");

  // All data hangs under the group's name
  if (!m_sidreRootGroup->hasGroup(name))
  {
    return nullptr;
  }
  axom::sidre::Group* group = m_sidreRootGroup->getGroup(name);
  if (group == nullptr)
  {
    return nullptr;
  }

  if (!group->hasView("value"))
  {
    return nullptr;
  }
  return group->getView("value");
}

bool Inlet::get(const std::string& name, bool& value)
{
  axom::sidre::View* valueView = baseGet(name);
  if (valueView == nullptr)
  {
    return false;
  }

  // There is no boolean type in conduit/sidre so we use int8
  if (valueView->getTypeID() != axom::sidre::INT8_ID)
  {
    std::string msg = fmt::format("Boolean named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  int8 intValue = valueView->getScalar();
  if (intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format("Invalid integer value stored in boolean"
                                  " value named {0}",
                                  name);
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  value = (bool)intValue;
  return true;
}

bool Inlet::get(const std::string& name, double& value)
{
  axom::sidre::View* valueView = baseGet(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::DOUBLE_ID)
  {
    std::string msg = fmt::format("Double named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  value = valueView->getScalar();
  return true;
}

bool Inlet::get(const std::string& name, int& value)
{
  axom::sidre::View* valueView = baseGet(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::INT_ID)
  {
    std::string msg = fmt::format("Integer named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  value = valueView->getScalar();
  return true;
}

bool Inlet::get(const std::string& name, std::string& value)
{
  axom::sidre::View* valueView = baseGet(name);
  if (valueView == nullptr)
  {
    return false;
  }

  if (valueView->getTypeID() != axom::sidre::CHAR8_STR_ID)
  {
    std::string msg = fmt::format("String named '{0}' was asked for"
                                  " but recieved type {1}",
                                  name, valueView->getTypeID());
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  const char* valueStr = valueView->getString();
  if (valueStr == nullptr)
  {
    value = std::string("");
  }
  value = std::string(valueStr);
  return true;
}

void Inlet::registerDocWriter(std::shared_ptr<DocWriter> writer) {
  m_docWriter = writer;
}

void Inlet::writeDoc() {
  if (m_docEnabled) {
    m_docWriter->writeDocumentation();
  }
}

bool Inlet::verify() {
  bool verifySuccess = true;
  verifyRecursive(m_sidreRootGroup, verifySuccess);
  return verifySuccess;
}

void Inlet::verifyRecursive(axom::sidre::Group* sidreGroup, bool& verifySuccess) {
  SLIC_ASSERT_MSG(sidreGroup, "Root was nullptr");
  if (sidreGroup == m_sidreRootGroup && sidreGroup->hasView("warningFlag")) {
    verifySuccess = false;
  }

  if (sidreGroup->hasView("required")) {
    int8 required = sidreGroup->getView("required")->getData();
    if (required && !sidreGroup->hasView("value")) {
      std::string msg = fmt::format("Inlet: {0}: Required field was not specified in Input Deck", 
                                    sidreGroup->getPathName());
      issueWarning(msg, m_sidreRootGroup);
      verifySuccess = false;
    }
  }

  if (sidreGroup->hasView("value")) {
    if (sidreGroup->hasView("range")) {
      auto type = sidreGroup->getView("range")->getTypeID();
      if (type == axom::sidre::INT_ID) {
        int* range = sidreGroup->getView("range")->getArray();
        int val = sidreGroup->getView("value")->getScalar();
        // Checks if the value is not within the specified range
        if (!(range[0] <= val && val <= range[1])) {
          verifySuccess = false;
        }
      } else if (type == axom::sidre::DOUBLE_ID) {
        double* range = sidreGroup->getView("range")->getArray();
        double val = sidreGroup->getView("value")->getScalar();
        if (!(range[0] <= val && val <= range[1])) {
          verifySuccess = false;
        }
      } else {
        verifySuccess = false;
      }
    } else if (sidreGroup->hasView("validValues")) {
      int val = sidreGroup->getView("value")->getScalar();
      int* range = sidreGroup->getView("validValues")->getArray();
      size_t size = sidreGroup->getView("validValues")->getBuffer()->getNumElements();
      // Checks if the value specified is one of the valid values
      bool found = false;
      for (size_t i = 0; i < size; i++) {
        if (range[i] == val) {
          found = true;
          break;
        }
      }
      if (!found) {
        verifySuccess = false;
      }
    }
  }
  
  axom::sidre::IndexType groupIndex = sidreGroup->getFirstValidGroupIndex();
  while (axom::sidre::indexIsValid(groupIndex)) {
    verifyRecursive(sidreGroup->getGroup(groupIndex), verifySuccess);
    groupIndex = sidreGroup->getNextValidGroupIndex(groupIndex);
  }
}

} // end namespace inlet
} // end namespace axom
