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

#include <memory>
#include <functional>
#include <vector>
#include <initializer_list>

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
class VerifiableScalar
{
public:
  /*!
   *****************************************************************************
   * \brief Set the required status of this object.
   *
   * Set whether this object is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether object is required
   *
   * \return Shared pointer to calling object, for chaining
   *****************************************************************************
   */
  virtual std::shared_ptr<VerifiableScalar> required(bool isRequired = true) = 0;

  /*!
   *****************************************************************************
   * \brief Return the required status.
   *
   * Return that this object is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this object is required
   *****************************************************************************
   */
  virtual bool isRequired() const = 0;

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this object's contents
   * during the verification stage.
   * 
   * \param [in] The function object.
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> registerVerifier(
    std::function<bool(const axom::inlet::Field&)> lambda) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> defaultValue(
    const std::string& value) = 0;
  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default string value
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> defaultValue(const char* value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default boolean value
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> defaultValue(bool value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default integer value
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> defaultValue(int value) = 0;

  /*!
   *****************************************************************************
   * \brief Set the default value of this object.
   *
   * Set the default value for the object in the input file.
   *
   * \param [in] value The default double value
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> defaultValue(double value) = 0;

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
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> range(double startVal,
                                                  double endVal) = 0;

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
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> range(int startVal, int endVal) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed integer values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::vector<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An vector containing the set of allowed double values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::vector<double>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set A vector containing the set of allowed string values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::vector<std::string>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the set of allowed C-string 
   * values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<const char*>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values for this object.
   *
   * \param [in] set An initializer list containing the valid integer values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<int>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Set the valid values.
   *
   * \param [in] set An initializer list containing the valid double values
   *
   * \return Shared pointer to calling object for chaining
   *****************************************************************************
  */
  virtual std::shared_ptr<VerifiableScalar> validValues(
    const std::initializer_list<double>& set) = 0;

  /*!
   *****************************************************************************
   * \brief Verifies the object to make sure it satisfies the imposed requirements
   *****************************************************************************
  */
  virtual bool verify() const = 0;
};

}  // namespace inlet
}  // namespace axom

#endif  // INLET_VERIFIABLE_SCALAR_HPP
