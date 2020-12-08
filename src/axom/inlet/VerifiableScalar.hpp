// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file VerifiableScalar.hpp
 *
 * \brief This file defines an interface for scalars that are verifiable
 *******************************************************************************
 */

#ifndef INLET_VERIFIABLE_SCALAR_HPP
#define INLET_VERIFIABLE_SCALAR_HPP

#include "axom/inlet/Verifiable.hpp"

namespace axom
{
namespace inlet
{
// Forward declarations
class Field;

/*!
 *******************************************************************************
 * \class VerifiableScalar
 *
 * \brief Basic interface for verifiable scalar values of Inlet primitive type,
 * namely int, double, bool, or std::string - implementations can use this
 * directly (inlet::Field) or forward to all elements of a container (inlet::AggregrateField)
 * 
 * In practice this interface is used for the Field and AggregateField types.
 * In addition to the ability to mark things as required and use user-defined lambdas
 * (as in Verifiable), default values can be provided, as can ranges of valid
 * values and discrete sets of valid values.
 *******************************************************************************
 */
class VerifiableScalar : public Verifiable<Field>
{
public:
  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(const std::string& value) = 0;
  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(const char* value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default boolean value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(bool value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default integer value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(int value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default double value
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& defaultValue(double value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the range of this object.
   *
   * Set the continuous range for the object in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& range(double startVal, double endVal) = 0;

  /*!
   *****************************************************************************
   * \brief Set the range of this object.
   *
   * Set the continuous range for the object in the input file.
   *
   * \param [in] startVal The start of the range
   * 
   * \param [in] endVal The end of the range
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& range(int startVal, int endVal) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<double>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::vector<std::string>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the set of allowed C-string 
   * values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(
    const std::initializer_list<const char*>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::initializer_list<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Reference to calling object for chaining
   *****************************************************************************
  */
  virtual VerifiableScalar& validValues(const std::initializer_list<double>& set) = 0;
};

}  // namespace inlet
}  // namespace axom

#endif  // INLET_VERIFIABLE_SCALAR_HPP
