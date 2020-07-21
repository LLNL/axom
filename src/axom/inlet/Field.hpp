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

namespace axom
{
namespace inlet
{

// enum class FieldType {BOOL, INT, DOUBLE, STRING, UNSPECIFIED};

/*!
 *******************************************************************************
 * \class Field
 *
 * \brief Provides functions to help define how individual field variables in an
 *        input deck are expected to behave.  It also holds the Sidre Group to 
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
   * \param [in] type FieldType specifying the data type of this Field instance.
   * Default is FieldType::UNSPECIFIED.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for Input Deck this Field instance belongs to.
   *****************************************************************************
   */
  Field(axom::sidre::Group* sidreGroup, 
        axom::sidre::DataTypeId type = axom::sidre::DataTypeId::NO_TYPE_ID,
        bool docEnabled = true) :
        m_sidreGroup(sidreGroup), m_type(type), m_docEnabled(docEnabled) {}

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
  std::shared_ptr<Field> addDefaultValue(const std::string& value);

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
  std::shared_ptr<Field> addDefaultValue(const char* value);

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
  std::shared_ptr<Field> addDefaultValue(bool value);

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
  std::shared_ptr<Field> addDefaultValue(int value);

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
  std::shared_ptr<Field> addDefaultValue(double value);

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
   * \brief Set the range of this Field.
   *
   * Set the continuous range for the Field in the input deck.
   *
   * \param [in] set An integer array indicating the set of allowed values
   * 
   * \param [in] size The number of elements in the array
   *
   * \return Shared pointer to this Field instance
   *****************************************************************************
  */
  std::shared_ptr<Field> validValues(int* set, size_t size);

private:
  // This Field's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
  axom::sidre::DataTypeId m_type = axom::sidre::DataTypeId::NO_TYPE_ID;
  bool m_docEnabled = false;
};

} // end namespace inlet
} // end namespace axom

#endif
