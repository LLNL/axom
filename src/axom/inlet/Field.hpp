// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Field.hpp
 *
 * \brief This file contains the class definition of Inlet's Field and AggregateField classes.
 *******************************************************************************
 */

#ifndef INLET_FIELD_HPP
#define INLET_FIELD_HPP

#include "axom/sidre.hpp"
#include "axom/inlet/Verifiable.hpp"

#include <memory>
#include <type_traits>
#include <functional>

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \enum InletType
 *
 * \brief Enumeration of basic types for things in inlet
 *******************************************************************************
 */
enum class InletType
{
  Nothing,
  Bool,
  String,
  Integer,
  // TODO: Unsigned integer
  Double,
  Object,
  Array
};

/*!
 *******************************************************************************
 * \class Field
 *
 * \brief Provides functions to help define how individual field variables in an
 * input file are expected to behave.  It also holds the Sidre Group to 
 * the individual field.
 *
 * \see Inlet Table
 *******************************************************************************
 */
class Field : public std::enable_shared_from_this<Field>, public VerifiableScalar
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Field class.
   *
   * This class provides functions to define the behavior of the Field
   * data already read and stored in the given Sidre Group.
   *
   * \param [in] sidreGroup Pointer to the already created Sidre Group.
   * \param [in] root Pointer to the sidreRootGroup containing this Field
   * \param [in] type FieldType specifying the data type of this Field instance.
   * Default is FieldType::UNSPECIFIED.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for Input file this Field instance belongs to.
   *****************************************************************************
   */
  Field(axom::sidre::Group* sidreGroup,
        axom::sidre::Group* root,
        axom::sidre::DataTypeId type = axom::sidre::DataTypeId::NO_TYPE_ID,
        bool docEnabled = true)
    : m_sidreGroup(sidreGroup)
    , m_sidreRootGroup(root)
    , m_type(type)
    , m_docEnabled(docEnabled)
  { }

  virtual ~Field() = default;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Field.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Field class.
   *
   * \return Pointer to the Sidre Group class for this Field
   *****************************************************************************
   */
  axom::sidre::Group* sidreGroup() { return m_sidreGroup; };

  /*!
   *****************************************************************************
   * \brief Set the required status of this Field.
   *
   * Set whether this Field is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Field is required
   *
   * \return Shared pointer to this instance of this class
   *****************************************************************************
   */
  std::shared_ptr<VerifiableScalar> required(bool isRequired);

  /*!
   *****************************************************************************
   * \brief Return the required status of this Field.
   *
   * Return that this Field is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this Field is required
   *****************************************************************************
   */
  bool isRequired();

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(const std::string& value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(const char* value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default boolean value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(bool value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default integer value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(int value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default double value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(double value);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> range(double startVal, double endVal);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> range(int startVal, int endVal);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<double>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<std::string>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the set of allowed C-string 
   * values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<const char*>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<double>& set);

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this Field's contents
   * during the verification stage.
   * 
   * \param [in] The function object that will be called by Field::verify().
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> registerVerifier(
    std::function<bool(Proxy&)> lambda);

  /*!
   *****************************************************************************
   * \brief Called by Inlet::verify to verify the contents of this Field.
   *****************************************************************************
  */
  bool verify();

  /*!
   *****************************************************************************
   * \return The full name for this Field.
   *****************************************************************************
  */
  std::string name();

  /*!
   *****************************************************************************
   * \brief Returns a value of primitive type
   * 
   * \return The value
   * \tparam T The type to retrieve
   *****************************************************************************
   */
  template <typename T>
  T get();

  /*!
   *****************************************************************************
   * \brief Returns the type of the stored value
   * 
   * \return The type
   * \see InletType
   *****************************************************************************
   */
  InletType type() const;

private:
  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set A vector containing the set of allowed scalar values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  template <typename T>
  void setScalarValidValues(std::vector<T> set);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  template <typename T>
  void setRange(T startVal, T endVal);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  template <typename T>
  void setDefaultValue(T value);

  /*!
   *****************************************************************************
   * \brief Checks the validity of a field value
   *
   * \param [in] view The view to verify
   *
   * \return Whether the value satisfied all constraints
   *****************************************************************************
  */
  bool verifyValue(axom::sidre::View& view);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is within the range.
   * 
   * \param [in] view The view containing the value that will be checked.
   * 
   * \return true if the given value was within its respective range, else false.
   * \pre T must define bool operator<=(T, T)
   *****************************************************************************
   */
  template <typename T>
  bool checkRange(axom::sidre::View& view);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is found in the list of valid values.
   * 
   * \param [in] view The view containing the value that will be checked.
   * 
   * \return true if the given target was found in its respective valid values, 
   *  else false.
   *****************************************************************************
   */
  template <typename T>
  bool searchValidValues(axom::sidre::View& view);

  /*!
   *****************************************************************************
   * \brief Checks the existence and type of the value for the field
   *
   * \param [in] expected The expected type for the value
   *
   * \return Non-owning pointer to the Sidre view containing the value
   * \note Treats a nonexistent value or type mismatch as an error and will
   * emit a SLIC_ERROR accordingly
   *****************************************************************************
  */
  axom::sidre::View* checkExistenceAndType(const axom::sidre::DataTypeId expected);

  // This Field's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  axom::sidre::DataTypeId m_type = axom::sidre::DataTypeId::NO_TYPE_ID;
  bool m_docEnabled;
  std::function<bool(Proxy&)> m_verifier;
};

// Prototypes for template specializations
template <>
bool Field::get<bool>();

template <>
int Field::get<int>();

template <>
double Field::get<double>();

template <>
std::string Field::get<std::string>();

template <>
inline bool Field::searchValidValues<std::string>(axom::sidre::View& view);

/*!
   *****************************************************************************
   * \brief A wrapper class that enables constraints on groups of Fields
   *****************************************************************************
  */
class AggregateField : public std::enable_shared_from_this<AggregateField>,
                       public VerifiableScalar
{
public:
  AggregateField(std::vector<std::shared_ptr<VerifiableScalar>>&& fields)
    : m_fields(std::move(fields))
  { }

  virtual ~AggregateField() = default;

  /*!
   *****************************************************************************
   * \brief Called by Inlet::verify to verify the contents of this Field.
   *****************************************************************************
  */
  bool verify();

  /*!
   *****************************************************************************
   * \brief Set the required status of this Field.
   *
   * Set whether this Field is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Field is required
   *
   * \return Shared pointer to this instance of this class
   *****************************************************************************
   */
  std::shared_ptr<VerifiableScalar> required(bool isRequired);

  /*!
   *****************************************************************************
   * \brief Return the required status of this Field.
   *
   * Return that this Field is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this Field is required
   *****************************************************************************
   */
  bool isRequired();

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(const std::string& value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(const char* value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default boolean value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(bool value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default integer value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(int value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input file.
   *
   * \param [in] value The default double value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> defaultValue(double value);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> range(double startVal, double endVal);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> range(int startVal, int endVal);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<double>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(const std::vector<std::string>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the set of allowed C-string 
   * values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<const char*>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<double>& set);
  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this Field's contents
   * during the verification stage.
   * 
   * \param [in] The function object that will be called by Field::verify().
   *****************************************************************************
  */
  std::shared_ptr<VerifiableScalar> registerVerifier(
    std::function<bool(Proxy&)> lambda);

private:
  std::vector<std::shared_ptr<VerifiableScalar>> m_fields;
};

}  // end namespace inlet
}  // end namespace axom

#endif
