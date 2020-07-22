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

std::shared_ptr<Field> Field::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr, "Inlet's Field specific Sidre Datastore Group not set");

  if (m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format("Inlet Field has already defined required value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
  if (m_type != axom::sidre::DataTypeId::INT_ID) {
    SLIC_WARNING("Field value type did not match INT");
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
  if (m_type != axom::sidre::DataTypeId::DOUBLE_ID) {
    SLIC_WARNING("Field value type did not match DOUBLE");
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }
  if (m_sidreGroup->hasView("defaultValue")) {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }

  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
  if (m_type != axom::sidre::DataTypeId::INT_ID) {
    SLIC_WARNING("Field value type did not match INT");
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }
  
  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
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
  if (m_type != axom::sidre::DataTypeId::INT_ID) {
    SLIC_WARNING("Field value type did not match INT");
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  }

  if (m_sidreGroup->hasView("range") || m_sidreGroup->hasView("validValues")) {
    std::string msg = fmt::format("Inlet Field has already defined range: {0}",
                                   m_sidreGroup->getName());
    SLIC_WARNING(msg);
    if (!m_sidreRootGroup->hasView("warningFlag")) {
      m_sidreRootGroup->createViewScalar("warningFlag", 1);
    }
  } else {
    auto view = m_sidreGroup->createViewAndAllocate("validValues", 
                                                    axom::sidre::INT_ID, set.size());
    view->getBuffer()->copyBytesIntoBuffer(&set[0], set.size()*sizeof(int));
  }
  return shared_from_this();
}



} // end namespace inlet
} // end namespace axom
