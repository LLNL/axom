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
#include "axom/inlet/inlet_utils.hpp"
#include <algorithm>

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
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);    
    return false;
  }

  int8 intValue = valueView->getScalar();
  if (intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format("Invalid integer value stored in boolean"
                                  " value named {0}",
                                  name);
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);    
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
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
      std::string msg = fmt::format("Inlet: {0}: Required field not specified", 
                                    sidreGroup->getPathName());
      SLIC_WARNING(msg);
      setWarningFlag(m_sidreRootGroup);
      verifySuccess = false;
    }
  }
  if (sidreGroup->hasView("value") && !verifyValue(sidreGroup)) {
    verifySuccess = false;
    std::string msg = fmt::format("Inlet: {0}: Value did not meet range/valid "
                                  "value(s) constraints", sidreGroup->getPathName());
    SLIC_WARNING(msg);
  }
  if (sidreGroup->hasView("defaultValue") && !verifyDefaultValue(sidreGroup)) {
    verifySuccess = false;
    std::string msg = fmt::format("Inlet: {0}: Default value did not meet range/valid "
                                  "value(s) constraints", sidreGroup->getPathName());
    SLIC_WARNING(msg);
  }
  
  axom::sidre::IndexType groupIndex = sidreGroup->getFirstValidGroupIndex();
  while (axom::sidre::indexIsValid(groupIndex)) {
    if (sidreGroup->getGroupName(groupIndex) != "validStringValues") {
      verifyRecursive(sidreGroup->getGroup(groupIndex), verifySuccess);
    }
    groupIndex = sidreGroup->getNextValidGroupIndex(groupIndex);
  }
}

bool Inlet::verifyValue(axom::sidre::Group* sidreGroup) {
  auto type = sidreGroup->getView("value")->getTypeID();
  if (sidreGroup->hasView("validValues")) {
    if (type == axom::sidre::INT_ID) {
      int val = sidreGroup->getView("value")->getScalar();
      return searchValidValues(sidreGroup, val);
    } else {
      double val = sidreGroup->getView("value")->getScalar();
      return searchValidValues(sidreGroup, val);
    }
  } else if (sidreGroup->hasView("range")) {
    if (type == axom::sidre::INT_ID) {
      int val = sidreGroup->getView("value")->getScalar();
      return checkRange(sidreGroup, val);
    } else {
      double val = sidreGroup->getView("value")->getScalar();
      return checkRange(sidreGroup, val);
    }
  } else if (sidreGroup->hasGroup("validStringValues")) {
    std::string val = sidreGroup->getView("value")->getString();
    return searchValidValues(sidreGroup->getGroup("validStringValues"), val);
  }
  return true;
}

bool Inlet::verifyDefaultValue(axom::sidre::Group* sidreGroup) {
  auto view = sidreGroup->getView("defaultValue");
  auto type = view->getTypeID();
  if (sidreGroup->hasView("validValues")) {
    if (type == axom::sidre::INT_ID) {
      int val = view->getScalar();
      return searchValidValues(sidreGroup, val);
    } else {
      double val = view->getScalar();
      return searchValidValues(sidreGroup, val);
    }
  } else if (sidreGroup->hasView("range")) {
    if (type == axom::sidre::INT_ID) {
      int val = view->getScalar();
      return checkRange(sidreGroup, val);
    } else {
      double val = view->getScalar();
      return checkRange(sidreGroup, val);
    }
  } else if (sidreGroup->hasGroup("validStringValues")) {
    std::string val = view->getString();
    return searchValidValues(sidreGroup->getGroup("validStringValues"), val);
  }
  return true;
}

bool Inlet::checkRange(axom::sidre::Group* sidreGroup, int val) {
  int* range = sidreGroup->getView("range")->getArray();
  return range[0] <= val && val <= range[1];
}

bool Inlet::checkRange(axom::sidre::Group* sidreGroup, double val) {
  double* range = sidreGroup->getView("range")->getArray();
  return range[0] <= val && val <= range[1];
}

bool Inlet::searchValidValues(axom::sidre::Group* sidreGroup, int target) {
  auto view = sidreGroup->getView("validValues");
  int* valuesArray = view->getArray();
  size_t size = view->getBuffer()->getNumElements();
  int* result = std::find(valuesArray, valuesArray + size, target);
  return result != valuesArray + size;
}

bool Inlet::searchValidValues(axom::sidre::Group* sidreGroup, double target) {
  auto view = sidreGroup->getView("validValues");
  double* valuesArray = view->getArray();
  size_t size = view->getBuffer()->getNumElements();
  double* result = std::find(valuesArray, valuesArray + size, target);
  return result != valuesArray + size;
}

bool Inlet::searchValidValues(axom::sidre::Group* sidreGroup, std::string value) {
  auto idx = sidreGroup->getFirstValidViewIndex();
  while(axom::sidre::indexIsValid(idx)) {
    if (sidreGroup->getView(idx)->getString() == value) {
      return true;
    }
    idx = sidreGroup->getNextValidViewIndex(idx);
  }
  return false;
}



} // end namespace inlet
} // end namespace axom
