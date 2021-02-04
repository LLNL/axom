// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Writer.hpp
 *
 * \brief This file contains the abstract base class definition of Writer.
 *******************************************************************************
 */

#ifndef INLET_WRITER_HPP
#define INLET_WRITER_HPP

namespace axom
{
namespace inlet
{
// Forward declaration
class Table;

/*!
 *******************************************************************************
 * \class Writer
 *
 * \brief Abstract base class defining the interface of all Writer
 *  classes.
 *
 *  Concrete instances need to inherit from this class and implement these
 *  functions.
 *
 * \see SphinxWriter
 *******************************************************************************
 */
class Writer
{
public:
  virtual ~Writer() = default;

  /*!
   *****************************************************************************
   * \brief Generates documentation for a Table and its child Fields/Functions
   * \param [in] table The Table to generate documentation for
   *
   * \note Implementers of this function are not responsible for generating
   * documentation for child Tables of this Table - only child Fields/Functions
   *****************************************************************************
   */
  virtual void documentTable(const Table& table) = 0;

  /*!
   *****************************************************************************
   * \brief Finalizes documentation generation (e.g., by writing it to a file)
   * 
   * This is a hint to implementers that no further tables will be documented
   *****************************************************************************
   */
  virtual void finalize() = 0;
};

}  // namespace inlet
}  // namespace axom

#endif
