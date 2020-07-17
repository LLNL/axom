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
    return false;
  }

  return (bool)intValue;
}

std::shared_ptr<Field> Field::addDefaultString(std::string value) {
  SLIC_ASSERT_MSG(m_type == FieldType::STRING, "Default value type did not match STRING");

  if (m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
  } else {
    m_sidreGroup->createViewString("defaultValue", value);
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::addDefaultBool(bool value) {
  SLIC_ASSERT_MSG(m_type == FieldType::BOOL, "Default value type did not match BOOL");
  if (m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
  } else {
    m_sidreGroup->createViewScalar("defaultValue", (int8)value);
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::addDefaultInt(int value) {
  SLIC_ASSERT_MSG(m_type == FieldType::INT, "Default value type did not match BOOL");
  if (m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
  } else {
    m_sidreGroup->createViewScalar("defaultValue", value);
  }
  return shared_from_this();
}

std::shared_ptr<Field> Field::addDefaultDouble(double value) {
  SLIC_ASSERT_MSG(m_type == FieldType::DOUBLE, "Default value type did not match DOUBLE");
  if (m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format("Inlet Field has already defined default value: {0}",
                                  m_sidreGroup->getName());
    SLIC_WARNING(msg);
  } else {
    m_sidreGroup->createViewScalar("defaultValue", value);
  }
  return shared_from_this();
}



} // end namespace inlet
} // end namespace axom
