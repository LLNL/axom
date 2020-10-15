// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
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

#include "axom/inlet/DocWriter.hpp"

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
   * \param [in] reader Shared pointer to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether documentation generation
   * is enabled. This also toggles the storing of documentation-specific information.
   *****************************************************************************
   */
  Inlet(std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true)
    : m_reader(reader)
    , m_sidreRootGroup(sidreRootGroup)
    , m_globalTable(
        std::make_shared<Table>("", "", m_reader, m_sidreRootGroup, docEnabled))
    , m_docEnabled(docEnabled)
  { }

  virtual ~Inlet() = default;

  /*!
   *****************************************************************************
   * \brief Returns the shared pointer to the Reader class.
   *
   * Provides access to the Reader class that is used to access the input file.
   *
   * \return Shared pointer to this instances' Reader class
   *****************************************************************************
   */
  std::shared_ptr<Reader> reader() { return m_reader; };

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
   * \brief Add a Table to the input file schema.
   *
   * Adds a Table to the input file schema. Tables hold a varying amount Fields
   * defined by the user.  By default, it is not required unless marked with
   * Table::isRequired(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the Table expected in the input file
   * \param [in] description Description of the Table
   *
   * \return Shared pointer to the created Table
   *****************************************************************************
   */
  std::shared_ptr<Table> addTable(const std::string& name,
                                  const std::string& description = "");

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
  std::shared_ptr<VerifiableScalar> addBool(const std::string& name,
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
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<VerifiableScalar> addDouble(
    const std::string& name,
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
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<VerifiableScalar> addInt(const std::string& name,
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
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<VerifiableScalar> addString(
    const std::string& name,
    const std::string& description = "");

  //
  // Functions that get the values out of the datastore
  //

  /*!
   *******************************************************************************
   * \brief Gets a value of arbitrary type out of the datastore
   * 
   * Retrieves a value of user-defined type, i.e., not double, int, bool, or string.
   * 
   * \param [in] name The name of the subtable representing the root of the object
   * \return The retrieved value
   * \tparam The type to retrieve
   * \pre Requires a specialization of FromInlet<T>
   * \note This function does not indicate failure in a way that can be handled
   * by a program - if an object of requested type does not exist at the specified
   * location, the program will terminate
   *******************************************************************************
   */
  template <typename T>
  T get(const std::string& name)
  {
    return m_globalTable->get<T>(name);
  }

  /*!
   *****************************************************************************
   * \brief Return whether a subobject with the given name is present in 
   * the datastore.
   *
   * \see Table::contains
   *****************************************************************************
   */
  bool contains(const std::string& name)
  {
    return m_globalTable->contains(name);
  }

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the datastore.
   * 
   * \see Table::operator[]
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) { return (*m_globalTable)[name]; }

  /*!
   *****************************************************************************
   * \brief Sets the associated DocWriter for the Inlet instance.
   *
   * Sets the associated DocWriter. If the DocWriter is already set, it will be
   * replaced by the one that was most recently set.
   *
   * \param [in] writer A pointer to a DocWriter object
   *
   *****************************************************************************
   */
  void registerDocWriter(std::shared_ptr<DocWriter> writer);

  /*!
   *****************************************************************************
   * \brief Writes input file documentation.
   *
   * This writes the input file's documentation through the registered DocWriter.
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
  bool verify();

  /*!
   *****************************************************************************
   * \return The global Table.
   *****************************************************************************
   */
  std::shared_ptr<Table> getGlobalTable() { return m_globalTable; }

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
  std::shared_ptr<Table> getTable(const std::string& name)
  {
    return m_globalTable->getTable(name);
  }

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Field.
   * 
   * \param [in] The string indicating the target name of the Field to be searched for.
   * 
   * \return The child Field matching the target name. If no such Field is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Field> getField(const std::string& name)
  {
    return m_globalTable->getField(name);
  }

  /*!
   *****************************************************************************
   * \brief Return whether a Table with the given name is present in Inlet.
   *
   * \return Boolean value indicating whether this Inlet contains the Table.
   *****************************************************************************
   */
  bool hasTable(const std::string& name)
  {
    return m_globalTable->hasTable(name);
  }

  /*!
   *****************************************************************************
   * \brief Return whether a Field with the given name is present in Inlet.
   *
   * \return Boolean value indicating whether this Inlet contains the Field.
   *****************************************************************************
   */
  bool hasField(const std::string& name)
  {
    return m_globalTable->hasField(name);
  }

  /*!
   *****************************************************************************
   * \return An unordered map from Field names to the child Field pointers for 
   * this Table.
   *****************************************************************************
   */
  std::unordered_map<std::string, std::shared_ptr<Field>> getChildFields()
  {
    return m_globalTable->getChildFields();
  }

  /*!
   *****************************************************************************
   * \return An unordered map from Table names to the child Table pointers for 
   * this Table.
   *****************************************************************************
   */
  std::unordered_map<std::string, std::shared_ptr<Table>> getChildTables()
  {
    return m_globalTable->getChildTables();
  }

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
  std::shared_ptr<Verifiable> addBoolArray(const std::string& name,
                                           const std::string& description = "")
  {
    return m_globalTable->addBoolArray(name, description);
  }

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
  std::shared_ptr<Verifiable> addIntArray(const std::string& name,
                                          const std::string& description = "")
  {
    return m_globalTable->addIntArray(name, description);
  }

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
  std::shared_ptr<Verifiable> addDoubleArray(const std::string& name,
                                             const std::string& description = "")
  {
    return m_globalTable->addDoubleArray(name, description);
  }

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
  std::shared_ptr<Verifiable> addStringArray(const std::string& name,
                                             const std::string& description = "")
  {
    return m_globalTable->addStringArray(name, description);
  }

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
  std::shared_ptr<Table> addGenericArray(const std::string& name,
                                         const std::string& description = "")
  {
    return m_globalTable->addGenericArray(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Get a boolean array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, bool>& map)
  {
    return m_globalTable->getArray(map);
  }

  /*!
   *****************************************************************************
   * \brief Get a int array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, int>& map)
  {
    return m_globalTable->getArray(map);
  }

  /*!
   *****************************************************************************
   * \brief Get a double array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, double>& map)
  {
    return m_globalTable->getArray(map);
  }

  /*!
   *****************************************************************************
   * \brief Get a string array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, std::string>& map)
  {
    return m_globalTable->getArray(map);
  }

  // TODO add update value functions
private:
  std::shared_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  std::shared_ptr<Table> m_globalTable;
  std::shared_ptr<DocWriter> m_docWriter;
  bool m_docEnabled;
};

}  // end namespace inlet
}  // end namespace axom

#endif
