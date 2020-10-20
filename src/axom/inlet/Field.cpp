// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/inlet/Field.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/Proxy.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
std::shared_ptr<VerifiableScalar> Field::required(bool isRequired)
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Field specific Sidre Datastore Group not set");

  if(m_sidreGroup->hasView("required"))
  {
    std::string msg = fmt::format(
      "[Inlet] Field has already defined "
      "required value: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return shared_from_this();
  }

  if(isRequired)
  {
    m_sidreGroup->createViewScalar("required", (int8)1);
  }
  else
  {
    m_sidreGroup->createViewScalar("required", (int8)0);
  }

  return shared_from_this();
}

bool Field::isRequired() const
{
  SLIC_ASSERT_MSG(m_sidreGroup != nullptr,
                  "[Inlet] Field specific Sidre Datastore Group not set");

  if(!m_sidreGroup->hasView("required"))
  {
    return false;
  }
  axom::sidre::View* valueView = m_sidreGroup->getView("required");
  if(valueView == nullptr)
  {
    return false;
  }
  int8 intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      "boolean value named {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
    return false;
  }

  return (bool)intValue;
}

template <typename T>
void Field::setDefaultValue(T value)
{
  if(m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format(
      "[Inlet] Field has already defined "
      "default value: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    if(m_docEnabled)
    {
      m_sidreGroup->createViewScalar("defaultValue", value);
    }
    if(!m_sidreGroup->hasView("value"))
    {
      m_sidreGroup->createViewScalar("value", value);
    }
  }
}

template <>
void Field::setDefaultValue<std::string>(std::string value)
{
  if(m_sidreGroup->hasView("defaultValue"))
  {
    std::string msg = fmt::format(
      "[Inlet] Field has already defined "
      "default value: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    if(m_docEnabled)
    {
      m_sidreGroup->createViewString("defaultValue", value);
    }
    if(!m_sidreGroup->hasView("value"))
    {
      m_sidreGroup->createViewString("value", value);
    }
  }
}

std::shared_ptr<VerifiableScalar> Field::defaultValue(const char* value)
{
  std::string str = "";
  if(value)
  {
    str = value;
  }
  setDefaultValue(std::string(value));
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::defaultValue(const std::string& value)
{
  if(m_type != axom::sidre::DataTypeId::CHAR8_STR_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue(value);
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::defaultValue(bool value)
{
  if(m_type != axom::sidre::DataTypeId::INT8_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match BOOL");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue((int8)value);
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::defaultValue(int value)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    setDefaultValue(value);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setDefaultValue((double)value);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT/DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::defaultValue(double value)
{
  if(m_type != axom::sidre::DataTypeId::DOUBLE_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
  }
  setDefaultValue(value);
  return shared_from_this();
}

template <typename T>
void Field::setRange(T startVal, T endVal)
{
  if(m_sidreGroup->hasView("range"))
  {
    std::string msg =
      fmt::format("[Inlet] Inlet Field has already defined range: {0}",
                  m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("validValues") ||
          m_sidreGroup->hasView("validStringValues"))
  {
    std::string msg = fmt::format(
      "[Inlet] Cannot set range for Inlet Field "
      "after setting valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto* view = m_sidreGroup->createViewAndAllocate("range", m_type, 2);
    T* pair = view->getArray();
    pair[0] = startVal;
    pair[1] = endVal;
  }
}

std::shared_ptr<VerifiableScalar> Field::range(int startVal, int endVal)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    this->setRange(startVal, endVal);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    this->setRange(static_cast<double>(startVal), static_cast<double>(endVal));
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT or DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::range(double startVal, double endVal)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::DOUBLE_ID:
    this->setRange(startVal, endVal);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return shared_from_this();
}

template <>
bool Field::get<bool>() const
{
  auto valueView = checkExistenceAndType(axom::sidre::INT8_ID);
  // There is no boolean type in conduit/sidre so we use int8
  int8 intValue = valueView->getScalar();
  if(intValue < 0 || intValue > 1)
  {
    std::string msg = fmt::format(
      "[Inlet] Invalid integer value stored in "
      " boolean value named {0}",
      name());
    SLIC_ERROR(msg);
  }

  return static_cast<bool>(intValue);
}

template <>
double Field::get<double>() const
{
  auto valueView = checkExistenceAndType(axom::sidre::DOUBLE_ID);
  return valueView->getScalar();
}

template <>
int Field::get<int>() const
{
  auto valueView = checkExistenceAndType(axom::sidre::INT_ID);
  return valueView->getScalar();
}

template <>
std::string Field::get<std::string>() const
{
  auto valueView = checkExistenceAndType(axom::sidre::CHAR8_STR_ID);

  const char* valueStr = valueView->getString();
  return (valueStr == nullptr) ? "" : valueStr;
}

axom::sidre::View* Field::checkExistenceAndType(
  const axom::sidre::DataTypeId expected) const
{
  axom::sidre::View* valueView = m_sidreGroup->getView("value");

  if(valueView == nullptr)
  {
    SLIC_ERROR("[Inlet] Field does not contain a value");
  }

  if(valueView->getTypeID() != expected)
  {
    std::string msg = fmt::format(
      "[Inlet] Field with name '{0}' was expected to be of type {1}"
      " but was actually of type {2}",
      name(),
      expected,
      valueView->getTypeID());
    SLIC_ERROR(msg);
  }

  return valueView;
}

InletType Field::type() const
{
  axom::sidre::View* valueView = m_sidreGroup->getView("value");
  if(valueView == nullptr)
  {
    return InletType::Nothing;
  }

  switch(valueView->getTypeID())
  {
  case axom::sidre::NO_TYPE_ID:
    return InletType::Nothing;
  case axom::sidre::INT8_ID:
    return InletType::Bool;
  case axom::sidre::INT_ID:
    return InletType::Integer;
  case axom::sidre::CHAR8_STR_ID:
    return InletType::String;
  case axom::sidre::DOUBLE_ID:
    return InletType::Double;
  default:
    std::string msg = fmt::format(
      "Type ID {0} for field not recognized, returning InletType::Nothing",
      valueView->getTypeID());
    SLIC_WARNING(msg);
    return InletType::Nothing;
  }
}

template <typename T>
void Field::setScalarValidValues(std::vector<T> set)
{
  if(m_sidreGroup->hasView("validValues") ||
     m_sidreGroup->hasView("validStringValues"))
  {
    std::string msg = fmt::format(
      "[Inlet] Inlet Field has already "
      "defined valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("range"))
  {
    std::string msg = fmt::format(
      "[Inlet] Cannot set valid values "
      "after defining range: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto view =
      m_sidreGroup->createViewAndAllocate("validValues", m_type, set.size());
    view->getBuffer()->copyBytesIntoBuffer(&set[0], set.size() * sizeof(T));
  }
}

std::shared_ptr<VerifiableScalar> Field::validValues(const std::vector<int>& set)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::INT_ID:
    setScalarValidValues(set);
    break;
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setScalarValidValues(std::vector<double>(set.begin(), set.end()));
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match INT OR DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::validValues(const std::vector<double>& set)
{
  switch(m_type)
  {
  case axom::sidre::DataTypeId::DOUBLE_ID:
    setScalarValidValues(set);
    break;
  default:  // incompatible type
    SLIC_WARNING("[Inlet] Field value type did not match DOUBLE");
    setWarningFlag(m_sidreRootGroup);
    break;
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::validValues(
  const std::vector<std::string>& set)
{
  if(m_type != axom::sidre::DataTypeId::CHAR8_STR_ID)
  {
    SLIC_WARNING("[Inlet] Field value type did not match STRING");
    setWarningFlag(m_sidreRootGroup);
  }
  if(m_sidreGroup->hasView("validValues") ||
     m_sidreGroup->hasView("validStringValues"))
  {
    std::string msg = fmt::format(
      "[Inlet] Inlet Field has already defined "
      "valid values: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else if(m_sidreGroup->hasView("range"))
  {
    std::string msg = fmt::format(
      "[Inlet] Cannot set valid values after defining "
      "range: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    setWarningFlag(m_sidreRootGroup);
  }
  else
  {
    auto group = m_sidreGroup->createGroup("validStringValues", true);
    for(std::string str : set)
    {
      group->createViewString("", str);
    }
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> Field::validValues(
  const std::initializer_list<const char*>& set)
{
  return validValues(std::vector<std::string>(set.begin(), set.end()));
}

std::shared_ptr<VerifiableScalar> Field::validValues(
  const std::initializer_list<int>& set)
{
  return validValues(std::vector<int>(set));
}

std::shared_ptr<VerifiableScalar> Field::validValues(
  const std::initializer_list<double>& set)
{
  return validValues(std::vector<double>(set));
}

std::shared_ptr<VerifiableScalar> Field::registerVerifier(
  std::function<bool(const Field&)> lambda)
{
  SLIC_WARNING_IF(m_verifier,
                  fmt::format("[Inlet] Verifier for Field "
                              "already set: {0}",
                              m_sidreGroup->getPathName()));
  m_verifier = lambda;
  return shared_from_this();
}

bool Field::verify() const
{
  // If this field was required, make sure soemething was defined in it
  if(m_sidreGroup->hasView("required"))
  {
    int8 required = m_sidreGroup->getView("required")->getData();
    if(required && !m_sidreGroup->hasView("value"))
    {
      std::string msg = fmt::format(
        "[Inlet] Required Field not "
        "specified: {0}",
        m_sidreGroup->getPathName());
      SLIC_WARNING(msg);
      return false;
    }
  }

  // Verify the provided value
  if(m_sidreGroup->hasView("value") &&
     !verifyValue(*m_sidreGroup->getView("value")))
  {
    std::string msg = fmt::format(
      "[Inlet] Value did not meet range/valid "
      "value(s) constraints: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    return false;
  }

  // Verify the default value
  if(m_sidreGroup->hasView("defaultValue") &&
     !verifyValue(*m_sidreGroup->getView("defaultValue")))
  {
    std::string msg = fmt::format(
      "[Inlet] Default value did not meet range/valid "
      "value(s) constraints: {0}",
      m_sidreGroup->getPathName());
    SLIC_WARNING(msg);
    return false;
  }

  // Lambda verification step
  if(m_verifier && !m_verifier(*this))
  {
    SLIC_WARNING(fmt::format("[Inlet] Field failed lambda verification: {0}",
                             m_sidreGroup->getPathName()));
    return false;
  }
  return true;
}

bool Field::verifyValue(axom::sidre::View& view) const
{
  auto type = view.getTypeID();
  if(m_sidreGroup->hasView("validValues"))
  {
    if(type == axom::sidre::INT_ID)
    {
      return searchValidValues<int>(view);
    }
    else
    {
      return searchValidValues<double>(view);
    }
  }
  else if(m_sidreGroup->hasView("range"))
  {
    if(type == axom::sidre::INT_ID)
    {
      return checkRange<int>(view);
    }
    else
    {
      return checkRange<double>(view);
    }
  }
  else if(m_sidreGroup->hasGroup("validStringValues"))
  {
    return searchValidValues<std::string>(view);
  }
  return true;
}

template <typename T>
bool Field::checkRange(axom::sidre::View& view) const
{
  T val = view.getScalar();
  T* range = m_sidreGroup->getView("range")->getArray();
  return range[0] <= val && val <= range[1];
}

template <typename T>
bool Field::searchValidValues(axom::sidre::View& view) const
{
  T target = view.getScalar();
  auto valid_vals = m_sidreGroup->getView("validValues");
  T* valuesArray = valid_vals->getArray();
  size_t size = valid_vals->getBuffer()->getNumElements();
  T* result = std::find(valuesArray, valuesArray + size, target);
  return result != valuesArray + size;
}

template <>
bool Field::searchValidValues<std::string>(axom::sidre::View& view) const
{
  auto string_group = m_sidreGroup->getGroup("validStringValues");
  std::string value = view.getString();
  auto idx = string_group->getFirstValidViewIndex();
  while(axom::sidre::indexIsValid(idx))
  {
    if(string_group->getView(idx)->getString() == value)
    {
      return true;
    }
    idx = string_group->getNextValidViewIndex(idx);
  }
  return false;
}

std::string Field::name() const
{
  return removePrefix(m_sidreRootGroup->getPathName(),
                      m_sidreGroup->getPathName());
}

bool AggregateField::verify() const
{
  return std::all_of(m_fields.begin(),
                     m_fields.end(),
                     [](const std::shared_ptr<VerifiableScalar>& field) {
                       return field->verify();
                     });
}

std::shared_ptr<VerifiableScalar> AggregateField::required(bool isRequired)
{
  for(auto& field : m_fields)
  {
    field->required(isRequired);
  }
  return shared_from_this();
}

bool AggregateField::isRequired() const
{
  for(auto& field : m_fields)
  {
    if(field->isRequired())
    {
      return true;
    }
  }
  return false;
}

std::shared_ptr<VerifiableScalar> AggregateField::defaultValue(
  const std::string& value)
{
  for(auto& field : m_fields)
  {
    field->defaultValue(value);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::defaultValue(const char* value)
{
  for(auto& field : m_fields)
  {
    field->defaultValue(value);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::defaultValue(bool value)
{
  for(auto& field : m_fields)
  {
    field->defaultValue(value);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::defaultValue(int value)
{
  for(auto& field : m_fields)
  {
    field->defaultValue(value);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::defaultValue(double value)
{
  for(auto& field : m_fields)
  {
    field->defaultValue(value);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::range(double startVal,
                                                        double endVal)
{
  for(auto& field : m_fields)
  {
    field->range(startVal, endVal);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::range(int startVal, int endVal)
{
  for(auto& field : m_fields)
  {
    field->range(startVal, endVal);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::vector<int>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::vector<double>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::vector<std::string>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::initializer_list<const char*>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::initializer_list<int>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::validValues(
  const std::initializer_list<double>& set)
{
  for(auto& field : m_fields)
  {
    field->validValues(set);
  }
  return shared_from_this();
}

std::shared_ptr<VerifiableScalar> AggregateField::registerVerifier(
  std::function<bool(const Field&)> lambda)
{
  for(auto& field : m_fields)
  {
    field->registerVerifier(lambda);
  }
  return shared_from_this();
}

}  // end namespace inlet
}  // end namespace axom
