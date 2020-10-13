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
#include <unordered_map>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Verifiable.hpp"

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
  virtual std::shared_ptr<Table> addTable(const std::string& name,
                                          const std::string& description) = 0;

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
  virtual std::shared_ptr<VerifiableScalar> addBool(
    const std::string& name,
    const std::string& description) = 0;

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
  virtual std::shared_ptr<VerifiableScalar> addDouble(
    const std::string& name,
    const std::string& description) = 0;

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
  virtual std::shared_ptr<VerifiableScalar> addInt(
    const std::string& name,
    const std::string& description) = 0;

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
  virtual std::shared_ptr<VerifiableScalar> addString(
    const std::string& name,
    const std::string& description) = 0;
  /*!
   *****************************************************************************
   * \brief Return whether a Table with the given name is present in this Table's subtree.
   *
   * \return Boolean value indicating whether this Table's subtree contains this Table.
   *****************************************************************************
   */
  virtual bool hasTable(const std::string& tableName) = 0;

  /*!
   *****************************************************************************
   * \brief Return whether a Field with the given name is present in this Table's
   *  subtree.
   *
   * \return Boolean value indicating whether this Table's subtree contains this Field.
   *****************************************************************************
   */
  virtual bool hasField(const std::string& fieldName) = 0;

  /*!
   *****************************************************************************
   * \return An unordered map from Field names to the child Field pointers for 
   * this Table.
   *****************************************************************************
   */
  virtual std::unordered_map<std::string, std::shared_ptr<Field>> getChildFields() = 0;

  /*!
   *****************************************************************************
   * \return An unordered map from Table names to the child Table pointers for 
   * this Table.
   *****************************************************************************
   */
  virtual std::unordered_map<std::string, std::shared_ptr<Table>> getChildTables() = 0;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Table.
   * 
   * \param [in] The string indicating the target name of the Table to be searched for.
   * 
   * \return The Table matching the target name. If no such Table is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  virtual std::shared_ptr<Table> getTable(const std::string& tableName) = 0;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Field.
   * 
   * \param [in] The string indicating the target name of the Field to be searched for.
   * 
   * \return The Field matching the target name. If no such Field is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  virtual std::shared_ptr<Field> getField(const std::string& fieldName) = 0;

  /*!
   *****************************************************************************
   * \brief Add an array of Boolean Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  virtual std::shared_ptr<Verifiable> addBoolArray(
    const std::string& name,
    const std::string& description = "",
    const std::string& path_override = "") = 0;

  /*!
   *****************************************************************************
   * \brief Add an array of Integer Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  virtual std::shared_ptr<Verifiable> addIntArray(
    const std::string& name,
    const std::string& description = "",
    const std::string& path_override = "") = 0;

  /*!
   *****************************************************************************
   * \brief Add an array of Double Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  virtual std::shared_ptr<Verifiable> addDoubleArray(
    const std::string& name,
    const std::string& description = "",
    const std::string& path_override = "") = 0;

  /*!
   *****************************************************************************
   * \brief Add an array of String Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  virtual std::shared_ptr<Verifiable> addStringArray(
    const std::string& name,
    const std::string& description = "",
    const std::string& path_override = "") = 0;

  /*!
   *****************************************************************************
   * \brief Add an array of Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  virtual std::shared_ptr<Table> addGenericArray(
    const std::string& name,
    const std::string& description = "") = 0;
};

}  // end namespace inlet
}  // end namespace axom

#endif
