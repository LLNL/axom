// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Verifiable.hpp
 *
 * \brief This file defines an interface for things that are verifiable
 *******************************************************************************
 */

#ifndef INLET_VERIFIABLE_HPP
#define INLET_VERIFIABLE_HPP

#include <functional>

namespace axom
{
namespace inlet
{
// Forward declarations
class Table;

/*!
 *******************************************************************************
 * \class Verifiable
 *
 * \brief Interface for trivially verifiable objects - namely those that can
 * be marked as required or checked with a user-provided lambda
 * 
 * In practice this interface is used for the Table and AggregateTable classes.
 * Currently the only supported means of verifying a composite type (table)
 * are the methods exposed by this interface.
 *******************************************************************************
 */
class Verifiable
{
public:
  // Should not be reassignable
  Verifiable& operator=(const Verifiable&) = delete;
  /*!
   *****************************************************************************
   * \brief Set the required status of this object.
   *
   * Set whether this object is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether object is required
   *
   * \return Reference to calling object, for chaining
   *****************************************************************************
   */
  virtual Verifiable& required(bool isRequired = true) = 0;

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
  virtual Verifiable& registerVerifier(
    std::function<bool(const axom::inlet::Table&)> lambda) = 0;

  /*!
   *****************************************************************************
   * \brief Verifies the object to make sure it satisfies the imposed requirements
   *****************************************************************************
  */
  virtual bool verify() const = 0;
};

}  // namespace inlet
}  // namespace axom

#endif  // INLET_VERIFIABLE_HPP
