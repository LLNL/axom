// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Field.hpp"
#include "axom/inlet/inlet_utils.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{

std::shared_ptr<Field> Field::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Field specific Sidre Datastore Group not set");

  if (m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format("Inlet Field has already defined required value: {0}",
                                  m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
                                  m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
  if (m_type != axom::sidre::DataTypeId::CHAR8_STR_ID) {
    SLIC_WARNING("Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
  if (m_type != axom::sidre::DataTypeId::INT8_ID) {
    SLIC_WARNING("Field value type did not match BOOL");
    setWarningFlag(m_sidreRootGroup);
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
  if (m_type != axom::sidre::DataTypeId::INT_ID 
      && m_type != axom::sidre::DataTypeId::DOUBLE_ID) {
    SLIC_WARNING("Field value type did not match INT");
    setWarningFlag(m_sidreRootGroup);
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  } else {
    if (m_docEnabled) {
      if (m_type == axom::sidre::DataTypeId::DOUBLE_ID) {
        m_sidreGroup->createViewScalar("defaultValue", (double)value);
      } else {
        m_sidreGroup->createViewScalar("defaultValue", value);
      }
    }
    if (!m_sidreGroup->hasView("value")) {
      if (m_type == axom::sidre::DataTypeId::DOUBLE_ID) {
        m_sidreGroup->createViewScalar("value", (double)value);
      } else {
        m_sidreGroup->createViewScalar("value", value);
      }
    }
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::defaultValue(double value) {
  if (m_type != axom::sidre::DataTypeId::DOUBLE_ID) {
    SLIC_WARNING("Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
  if (m_type != axom::sidre::DataTypeId::DOUBLE_ID) {
    SLIC_WARNING("Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
  }

  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
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
  if (m_type != axom::sidre::DataTypeId::INT_ID && m_type != axom::sidre::DataTypeId::DOUBLE_ID) {
    SLIC_WARNING("Field value type did not match INT");
    setWarningFlag(m_sidreRootGroup);
  }
  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")
      || m_sidreGroup->hasView("validStringValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  } else {
    if (m_type == axom::sidre::DataTypeId::INT_ID) {
      int tuple[2];
      tuple[0] = startVal;
      tuple[1] = endVal;
      auto view = m_sidreGroup->createViewAndAllocate("range", axom::sidre::INT_ID, 2);
      view->getBuffer()->copyBytesIntoBuffer(tuple, 2*sizeof(int));
    } else {
      double tuple[2];
      tuple[0] = startVal;
      tuple[1] = endVal;
      auto view = m_sidreGroup->createViewAndAllocate("range", axom::sidre::DOUBLE_ID, 2);
      view->getBuffer()->copyBytesIntoBuffer(tuple, 2*sizeof(double));
    }
    
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::validValues(std::vector<int> set) {
  if (m_type != axom::sidre::DataTypeId::INT_ID) {
    SLIC_WARNING("Field value type did not match INT");
    setWarningFlag(m_sidreRootGroup);
  }

  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")
      || m_sidreGroup->hasView("validStringValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  } else {
    auto view = m_sidreGroup->createViewAndAllocate("validValues", 
                                                    axom::sidre::INT_ID, set.size());
    view->getBuffer()->copyBytesIntoBuffer(&set[0], set.size()*sizeof(int));
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::validValues(const std::vector<std::string>& set) {
  if (m_type != axom::sidre::DataTypeId::CHAR8_STR_ID) {
    SLIC_WARNING("Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }

  if (m_sidreGroup->hasView("validStringValues") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined valid values: {0}",
                                   m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  } else {
    auto group = m_sidreGroup->createGroup("validStringValues", true);
    for (std::string str : set) {
      group->createViewString("", str);
    }
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::validValues(const std::initializer_list<const char*>& set) {
  return validValues(std::vector<std::string>(set.begin(), set.end()));
}

} // end namespace inlet
} // end namespace axom
