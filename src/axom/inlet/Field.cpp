// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Field.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

void issueWarning(bool condition, std::string msg, axom::sidre::Group* root) {
  if (!root) {
    SLIC_WARNING("Given sidre root was nullptr");
  } else if (condition) {
    SLIC_WARNING(msg);
    if (!root->hasView("warningFlag")) {
      root->createViewScalar("warningFlag", 1);
    }
  }
}

void issueWarning(std::string msg, axom::sidre::Group* root) {
  if (!root) {
    SLIC_WARNING("Given sidre root was nullptr");
  } else {
    SLIC_WARNING(msg);
    if (!root->hasView("warningFlag")) {
      root->createViewScalar("warningFlag", 1);
    }
  }
}

std::shared_ptr<Field> Field::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Field specific Sidre Datastore Group not set");

  if (m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format("Inlet Field has already defined required value: {0}",
                                  m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
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

bool Field::required()
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Field specific Sidre Datastore Group not set");

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
    issueWarning(msg, m_sidreRootGroup);
    return false;
  }

  return (bool)intValue;
}

std::shared_ptr<Field> Field::defaultValue(const char* value) {
  std::string str = "";
  if (value) {
    str = value;
  }
  return defaultValue(str);
}

std::shared_ptr<Field> Field::defaultValue(const std::string& value) {
  issueWarning(m_type != axom::sidre::DataTypeId::CHAR8_STR_ID, 
                    "Default value type did not match STRING", m_sidreRootGroup);
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    if (m_docEnabled) {
      m_sidreGroup->createViewString("defaultValue", value);
    }
    if (!m_sidreGroup->hasView("value")) {
      m_sidreGroup->createViewString("value", value);
    }
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::defaultValue(bool value) {
  issueWarning(m_type != axom::sidre::DataTypeId::INT8_ID, 
                    "Default value type did not match BOOL", m_sidreRootGroup);
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    if (m_docEnabled) {
      m_sidreGroup->createViewScalar("defaultValue", (int8)value);
    }
    if (!m_sidreGroup->hasView("value")) {
      m_sidreGroup->createViewScalar("value", (int8)value);
    }
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::defaultValue(int value) {
  issueWarning(m_type !=  axom::sidre::DataTypeId::INT_ID, 
                    "Default value type did not match INT", m_sidreRootGroup);
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    if (m_docEnabled) {
      m_sidreGroup->createViewScalar("defaultValue", value);
    }
    if (!m_sidreGroup->hasView("value")) {
      m_sidreGroup->createViewScalar("value", value);
    }
  }
  
  return shared_from_this();
}

std::shared_ptr<Field> Field::defaultValue(double value) {
  issueWarning(m_type != axom::sidre::DataTypeId::DOUBLE_ID, 
                    "Default value type did not match DOUBLE", m_sidreRootGroup);
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    if (m_docEnabled) {
      m_sidreGroup->createViewScalar("defaultValue", value);
    }
    if (!m_sidreGroup->hasView("value")) {
      m_sidreGroup->createViewScalar("value", value);
    }
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::range(double startVal, double endVal) {
  issueWarning(m_type != axom::sidre::DataTypeId::DOUBLE_ID, 
                       "Range value type did not match DOUBLE", m_sidreRootGroup);
  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    double tuple[2];
    tuple[0] = startVal;
    tuple[1] = endVal;
    auto view = m_sidreGroup->createViewAndAllocate("range", axom::sidre::DOUBLE_ID, 2);
    view->getBuffer()->copyBytesIntoBuffer(tuple, 2*sizeof(double));
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::range(int startVal, int endVal) {
  issueWarning(m_type !=  axom::sidre::DataTypeId::INT_ID, 
                       "Range value type did not match INT", m_sidreRootGroup);
  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    int tuple[2];
    tuple[0] = startVal;
    tuple[1] = endVal;
    auto view = m_sidreGroup->createViewAndAllocate("range", axom::sidre::INT_ID, 2);
    view->getBuffer()->copyBytesIntoBuffer(tuple, 2*sizeof(int));
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::validValues(std::vector<int> set) {
  issueWarning(m_type !=  axom::sidre::DataTypeId::INT_ID, 
               "Range value type did not match INT", m_sidreRootGroup);

  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    issueWarning(msg, m_sidreRootGroup);
  } else {
    auto view = m_sidreGroup->createViewAndAllocate("validValues", 
                                                    axom::sidre::INT_ID, set.size());
    view->getBuffer()->copyBytesIntoBuffer(&set[0], set.size()*sizeof(int));
  }
  return shared_from_this();
}



} // end namespace inlet
} // end namespace axom
