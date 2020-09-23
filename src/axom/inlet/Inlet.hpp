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

#include "axom/inlet/SchemaCreator.hpp"
#include "axom/inlet/Table.hpp"
#include "axom/inlet/Field.hpp"
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
class Inlet : public SchemaCreator
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
  std::shared_ptr<Field> addBool(const std::string& name,
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
  std::shared_ptr<Field> addDouble(const std::string& name,
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
  std::shared_ptr<Field> addInt(const std::string& name,
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
  std::shared_ptr<Field> addString(const std::string& name,
                                   const std::string& description = "");

  //
  // Functions that get the values out of the datastore
  //

  /*!
   *****************************************************************************
   * \brief Gets a Boolean value out of the Datastore.
   *
   * Retrieves the Field value out of the DataStore.  This Field may not have
   * been actually present in the input file and will be indicted by the return
   * value. 
   *
   * \param [in] name Name of the Field value to be gotten
   * \param [out] value Value to be filled
   *
   * \return True if the value was found in the Datastore
   *****************************************************************************
   */
  bool get(const std::string& name, bool& value);

  /*!
   *****************************************************************************
   * \brief Gets a Double value out of the Datastore.
   *
   * Retrieves the Field value out of the DataStore.  This Field may not have
   * been actually present in the input file and will be indicted by the return
   * value. 
   *
   * \param [in] name Name of the Field value to be gotten
   * \param [out] value Value to be filled
   *
   * \return True if the value was found in the Datastore
   *****************************************************************************
   */
  bool get(const std::string& name, double& value);

  /*!
   *****************************************************************************
   * \brief Gets a Integer value out of the Datastore.
   *
   * Retrieves the Field value out of the DataStore.  This Field may not have
   * been actually present in the input file and will be indicted by the return
   * value. 
   *
   * \param [in] name Name of the Field value to be gotten
   * \param [out] value Value to be filled
   *
   * \return True if the value was found in the Datastore
   *****************************************************************************
   */
  bool get(const std::string& name, int& value);

  /*!
   *****************************************************************************
   * \brief Gets a String value out of the Datastore.
   *
   * Retrieves the Field value out of the DataStore.  This Field may not have
   * been actually present in the input file and will be indicted by the return
   * value. 
   *
   * \param [in] name Name of the Field value to be gotten
   * \param [out] value Value to be filled
   *
   * \return True if the value was found in the Datastore
   *****************************************************************************
   */
  bool get(const std::string& name, std::string& value);

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
  std::unordered_map<std::string, std::shared_ptr<Field>> getChildFields() {
    return m_globalTable->getChildFields();
  }

  /*!
   *****************************************************************************
   * \return An unordered map from Table names to the child Table pointers for 
   * this Table.
   *****************************************************************************
   */
  std::unordered_map<std::string, std::shared_ptr<Table>> getChildTables() {
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
  std::shared_ptr<Table> addBoolArray(const std::string& name,
                                      const std::string& description = "") {
    return m_globalTable->addBoolArray(name,description);                                        
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
  std::shared_ptr<Table> addIntArray(const std::string& name,
                                     const std::string& description = "") {
    return m_globalTable->addIntArray(name,description);  
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
  std::shared_ptr<Table> addDoubleArray(const std::string& name,
                                        const std::string& description = "") {
    return m_globalTable->addDoubleArray(name,description);                                     
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
  std::shared_ptr<Table> addStringArray(const std::string& name,
                                        const std::string& description = "") {
    return m_globalTable->addStringArray(name,description);                                    
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
  bool getBoolArray(std::unordered_map<int, bool>& map) {
    return m_globalTable->getBoolArray(map);
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
  bool getIntArray(std::unordered_map<int, int>& map) {
    return m_globalTable->getIntArray(map);
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
  bool getDoubleArray(std::unordered_map<int, double>& map) {
    return m_globalTable->getDoubleArray(map);
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
  bool getStringArray(std::unordered_map<int, std::string>& map) {
    return m_globalTable->getStringArray(map);
  }

  // TODO add update value functions
private:
  /*!
   *****************************************************************************
   * \brief Verifies the contents of the sidreGroup according to Inlet 
   * requirements.
   *
   * This is the recursive internal helper for verify().
   * 
   * \param [in] sidreGroup The root of the sub-group to be verified.
   *
   * \param [out] verifySuccess Indicates whether the verification was 
   * successful: true if successful and false if not.
   *****************************************************************************
   */
  void verifyRecursive(axom::sidre::Group* sidreGroup, bool& verifySuccess);

  /*!
   *****************************************************************************
   * \brief Verifies the value of a Field.
   *
   * This checks whether the value in the Field is within the specified range or
   * contained in allowed values.
   * 
   * \param [in] sidreGroup The Sidre Group containing the value to be verified.
   *
   * \return boolean value indicating whether the verification was 
   * successful: true if successful and false if not.
   *****************************************************************************
   */
  bool verifyValue(axom::sidre::Group* sidreGroup);

  /*!
   *****************************************************************************
   * \brief Verifies the default value of a Field.
   *
   * This checks whether the default value in the Field is within the specified 
   * range or is one of the allowed values.
   * 
   * \param [in] sidreGroup The Sidre Group containing the default value to be 
   * verified.
   *
   * \return boolean value indicating whether the verification was 
   * successful: true if successful and false if not.
   *****************************************************************************
   */
  bool verifyDefaultValue(axom::sidre::Group* sidreGroup);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is within the range.
   * 
   * \param [in] sidreGroup The Sidre Group containing the range.
   * \param [in] value The integer value that will be checked.
   * 
   * \return true if the given value was within its respective range, else false.
   *****************************************************************************
   */
  bool checkRange(axom::sidre::Group* sidreGroup, int value);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is within the range.
   * 
   * \param [in] sidreGroup The Sidre Group containing the range.
   * \param [in] value The double value that will be checked.
   * 
   * \return true if the given value was within its respective range, else false.
   *****************************************************************************
   */
  bool checkRange(axom::sidre::Group* sidreGroup, double value);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is found in the list of valid values.
   * 
   * \param [in] sidreGroup The Sidre Group containing the valid values.
   * \param [in] value The target integer value that will be searched for.
   * 
   * \return true if the given target was found in its respective valid values, 
   *  else false.
   *****************************************************************************
   */
  bool searchValidValues(axom::sidre::Group* sidreGroup, int value);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is found in the list of valid values.
   * 
   * \param [in] sidreGroup The Sidre Group containing the valid values.
   * \param [in] value The target double value that will be searched for.
   * 
   * \return true if the given target was found in its respective valid values, 
   *  else false.
   *****************************************************************************
   */
  bool searchValidValues(axom::sidre::Group* sidreGroup, double value);

  /*!
   *****************************************************************************
   * \brief Checks if the given value is found in the list of valid values.
   * 
   * \param [in] sidreGroup The Sidre Group containing the valid values.
   * \param [in] value The target string value that will be searched for.
   * 
   * \return true if the given target was found in its respective list of valid
   * values, else false.
   *****************************************************************************
   */
  bool searchValidValues(axom::sidre::Group* sidreGroup, std::string value);

  axom::sidre::View* baseGet(const std::string& name);

  std::shared_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  std::shared_ptr<Table> m_globalTable;
  std::shared_ptr<DocWriter> m_docWriter;
  bool m_docEnabled;
};

}  // end namespace inlet
}  // end namespace axom

#endif
