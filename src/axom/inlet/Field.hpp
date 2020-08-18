// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Field.hpp
 *
 * \brief This file contains the class definition of Inlet's Field class.
 *******************************************************************************
 */

#ifndef INLET_FIELD_HPP
#define INLET_FIELD_HPP

#include "axom/sidre.hpp"

#include <memory>
#include <type_traits>

namespace axom
{
namespace inlet
{

/*!
 *******************************************************************************
 * \class Field
 *
 * \brief Provides functions to help define how individual field variables in an
 * input deck are expected to behave.  It also holds the Sidre Group to 
 * the individual field.
 *
 * \see Inlet Table
 *******************************************************************************
 */
class Field : public std::enable_shared_from_this<Field>
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
   * generation is enabled for Input Deck this Field instance belongs to.
   *****************************************************************************
   */
  Field(axom::sidre::Group* sidreGroup, axom::sidre::Group* root,
        axom::sidre::DataTypeId type = axom::sidre::DataTypeId::NO_TYPE_ID,
        bool docEnabled = true) : m_sidreGroup(sidreGroup), m_sidreRootGroup(root),
        m_type(type), m_docEnabled(docEnabled) {}

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
   * Set whether this Field is required, or not, to be in the input deck.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Field is required
   *
   * \return Shared pointer to this instance of this class
   *****************************************************************************
   */
  std::shared_ptr<Field> required(bool isRequired);

  /*!
   *****************************************************************************
   * \brief Return the required status of this Field.
   *
   * Return that this Field is required, or not, to be in the input deck.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this Field is required
   *****************************************************************************
   */
  bool required();

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> defaultValue(const std::string& value);

   /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> defaultValue(const char* value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default boolean value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> defaultValue(bool value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default integer value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> defaultValue(int value);

  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default double value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> defaultValue(double value);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input deck.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> range(double startVal, double endVal);

  /*!
   *****************************************************************************
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input deck.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> range(int startVal, int endVal);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(const std::vector<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(const std::vector<double>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(const std::vector<std::string>& set);

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
  std::shared_ptr<Field> validValues(const std::initializer_list<const char*>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(const std::initializer_list<int>& set);

  /*!
   *****************************************************************************
   * \brief Set the valid values for this Field.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(const std::initializer_list<double>& set);
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
   * Set the continuous range for the Field in the input deck.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  template<typename T>
  void setRange(T startVal, T endVal);
  
  /*!
   *****************************************************************************
   * \brief Set the default value of this Field.
   *
   * Set the default value for the Field in the input deck.
   *
   * \param [in] value The default value
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  template<typename T>
  void setDefaultValue(T value);

  // This Field's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  axom::sidre::DataTypeId m_type = axom::sidre::DataTypeId::NO_TYPE_ID;
  bool m_docEnabled = false;
};

} // end namespace inlet
} // end namespace axom

#endif
