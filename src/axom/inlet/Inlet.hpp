// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Inlet.hpp
 *
 * \brief This file contains the class definition of Inlet, the main class
 *        for the Inlet component.
 *******************************************************************************
 */

#ifndef INLET_INLET_HPP
#define INLET_INLET_HPP

#include <memory>
#include <string>
#include <vector>
#include <functional>

#include "axom/inlet/Table.hpp"
#include "axom/inlet/Field.hpp"
#include "axom/inlet/Proxy.hpp"
#include "axom/inlet/Reader.hpp"

#include "axom/sidre.hpp"

#include "axom/inlet/Writer.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class Inlet
 *
 * \brief This class is the main access point for all Inlet operations from
 *        from defining the schema of the users input file to getting the values
 *        out of the Sidre DataStore.
 *
 * \see Table Field
 *******************************************************************************
 */
class Inlet
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Inlet class.
   *
   * Creates an Inlet class that can then be used with the given Reader and will
   * store data under the given Sidre Group.
   *
   * \param [in] reader Unique (owning) pointer to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether documentation generation
   * is enabled. This also toggles the storing of documentation-specific information.
   *****************************************************************************
   */
  Inlet(std::unique_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true)
    : m_reader(std::move(reader))
    , m_sidreRootGroup(sidreRootGroup)
    , m_globalTable("", "", *m_reader, m_sidreRootGroup, docEnabled)
    , m_docEnabled(docEnabled)
  { }

  // Inlet objects must be move only - delete the implicit shallow copy constructor
  Inlet(const Inlet&) = delete;
  Inlet(Inlet&&) = default;

  virtual ~Inlet() = default;

  /*!
   *****************************************************************************
   * \brief Returns the reference to the Reader class.
   *
   * Provides access to the Reader class that is used to access the input file.
   *
   * \return Reference to this instances' Reader class
   *****************************************************************************
   */
  Reader& reader() { return *m_reader; };

  /*!
   *****************************************************************************
   * \brief Returns pointer to the root Sidre Group class for all of Inlet.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for all of Inlet.
   *
   * \return Pointer to the root Sidre Group for Inlet
   *****************************************************************************
   */
  axom::sidre::Group* sidreGroup() { return m_sidreRootGroup; };

  //
  // Functions that define the input file schema
  //

  /*!
   *****************************************************************************
   * \brief Add a structure to the input file schema.
   *
   * Adds a structure/record to the input file schema. Structures can contain
   * fields and/or substructures.  By default, it is not required unless marked with
   * Table::isRequired(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the struct expected in the input file
   * \param [in] description Description of the struct
   *
   * \return Reference to the created struct, as a Table
   *****************************************************************************
   */
  Table& addStruct(const std::string& name, const std::string& description = "");

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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addBool(const std::string& name,
                            const std::string& description = "");

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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addDouble(const std::string& name,
                              const std::string& description = "");

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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addInt(const std::string& name,
                           const std::string& description = "");

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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addString(const std::string& name,
                              const std::string& description = "");

  //
  // Functions that get the values out of the datastore
  //

  /*!
   *******************************************************************************
   * \brief Gets a value of arbitrary type out of the datastore
   * 
   * Retrieves a value of primitive or user-defined type.
   * 
   * \param [in] name The name of the subtable representing the root of the object
   * \return The retrieved value
   * \tparam The type to retrieve
   * \pre Requires a specialization of FromInlet<T> for user-defined types
   * \note This function does not indicate failure in a way that can be handled
   * by a program - if an object of requested type does not exist at the specified
   * location, the program will terminate
   *******************************************************************************
   */
  template <typename T>
  T get(const std::string& name) const
  {
    return m_globalTable.get<T>(name);
  }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name is present in 
   * the datastore.
   *
   * \see Table::contains
   *****************************************************************************
   */
  bool contains(const std::string& name) const
  {
    return m_globalTable.contains(name);
  }

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the datastore.
   * 
   * \see Table::operator[]
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const
  {
    return m_globalTable[name];
  }

  /*!
   *****************************************************************************
   * \brief Sets the associated Writer for the Inlet instance.
   *
   * Sets the associated Writer. If the Writer is already set, it will be
   * replaced by the one that was most recently set.
   *
   * \param [in] writer An owning pointer to a Writer object
   *
   *****************************************************************************
   */
  void registerWriter(std::unique_ptr<Writer> writer);

  /*!
   *****************************************************************************
   * \brief Writes input file documentation.
   *
   * This writes the input file's documentation through the registered Writer.
   *
   *****************************************************************************
   */
  void writeDoc();

  /*!
   *****************************************************************************
   * \brief Verifies the contents of the sidreGroup according to Inlet 
   * requirements.
   *
   * This recursively checks the correctness of each Field and Table in the Sidre
   * Group: ensuring that required Fields are specified, each Field's value 
   * and default value are within the specified range or are equal to a valid 
   * value, and types are consistent. Also ensures that the registered verification
   * functions hold true.
   * 
   * \return true if contents are correct and false if not.
   *
   *****************************************************************************
   */
  bool verify() const;

  /*!
   *****************************************************************************
   * \return The global Table.
   *****************************************************************************
   */
  Table& getGlobalTable() { return m_globalTable; }

  /*!
   *****************************************************************************
   * \brief Add an array of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addBoolArray(const std::string& name,
                                  const std::string& description = "")
  {
    return m_globalTable.addBoolArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addIntArray(const std::string& name,
                                 const std::string& description = "")
  {
    return m_globalTable.addIntArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of Double Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addDoubleArray(const std::string& name,
                                    const std::string& description = "")
  {
    return m_globalTable.addDoubleArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of String Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addStringArray(const std::string& name,
                                    const std::string& description = "")
  {
    return m_globalTable.addStringArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an array of user-defined type to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the array
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Table& addStructArray(const std::string& name,
                        const std::string& description = "")
  {
    return m_globalTable.addStructArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Get a function from the input deck
   *
   * \param [in] name        Name of the function
   * \param [in] ret_type    The return type of the function
   * \param [in] arg_types   The argument types of the function
   * \param [in] description Description of the Field
   *
   * \return Reference to the created Function
   *****************************************************************************
   */
  Verifiable<Function>& addFunction(const std::string& name,
                                    const FunctionTag ret_type,
                                    const std::vector<FunctionTag>& arg_types,
                                    const std::string& description = "")
  {
    return m_globalTable.addFunction(name, ret_type, arg_types, description);
  }
  /*!
   *****************************************************************************
   * \brief Add a dictionary of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Table>& addBoolDictionary(const std::string& name,
                                       const std::string& description = "")
  {
    return m_globalTable.addBoolDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Table>& addIntDictionary(const std::string& name,
                                      const std::string& description = "")
  {
    return m_globalTable.addIntDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of Double Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Table>& addDoubleDictionary(const std::string& name,
                                         const std::string& description = "")
  {
    return m_globalTable.addDoubleDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add a dictionary of String Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Table>& addStringDictionary(const std::string& name,
                                         const std::string& description = "")
  {
    return m_globalTable.addStringDictionary(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Add an dictionary of user-defined type to the input file schema.
   *
   * \param [in] name Name of the dictionary
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Table& addStructDictionary(const std::string& name,
                             const std::string& description = "")
  {
    return m_globalTable.addStructDictionary(name, description);
  }

  // TODO add update value functions
private:
  std::unique_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  Table m_globalTable;
  std::unique_ptr<Writer> m_writer;
  bool m_docEnabled;
};

}  // end namespace inlet
}  // end namespace axom

#endif
