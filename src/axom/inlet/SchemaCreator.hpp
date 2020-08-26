// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SchemaCreator.hpp
 *
 * \brief This file contains the pure virtual base class definition of
 *        SchemaCreator.
 *******************************************************************************
 */

#ifndef INLET_SCHEMACREATOR_HPP
#define INLET_SCHEMACREATOR_HPP

#include <memory>
#include <string>

#include "axom/inlet/Field.hpp"

namespace axom
{
namespace inlet
{

class Table;

/*!
 *******************************************************************************
 * \class SchemaCreator
 *
 * \brief Abstract base class defining the interface of all SchemaCreator
 *        classes.
 *
 *  Concrete instances need to inherit from this class and implement these
 *  functions. This ensures that Inlet and Table follow the same function
 *  signatures.
 *
 * \see Inlet Table
 *******************************************************************************
 */
class SchemaCreator
{
public:
  /*!
   *****************************************************************************
   * \brief Add a Table to the input file schema.
   *
   * Adds a Table to the input file schema. Tables hold a varying amount Fields
   * defined by the user.  By default, it is not required unless marked with
   * Table::required(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the Table expected in the input file
   * \param [in] description Description of the Table
   *
   * \return Shared pointer to the created Table
   *****************************************************************************
   */
  std::shared_ptr<Table> addTable(const std::string& name,
                                  const std::string& description);

  /*!
   *****************************************************************************
   * \brief Add a Boolean Field to the input file schema.
   *
   * Adds a Boolean Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addBool(const std::string& name,
                                 const std::string& description);

  /*!
   *****************************************************************************
   * \brief Add a Double Field to the input file schema.
   *
   * Adds a Double Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addDouble(const std::string& name,
                                   const std::string& description);

  /*!
   *****************************************************************************
   * \brief Add a Integer Field to the input file schema.
   *
   * Adds a Integer Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input file
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addInt(const std::string& name,
                                const std::string& description);

  /*!
   *****************************************************************************
   * \brief Add a String Field to the input file schema.
   *
   * Adds a String Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Table expected in the input file
   * \param [in] description Description of the Table
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addString(const std::string& name,
                                   const std::string& description);
};

} // end namespace inlet
} // end namespace axom

#endif
