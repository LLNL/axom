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
 *        from defining the schema of the users input deck to getting the values
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
   * \param [in] reader Shared pointer to the input deck Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   *****************************************************************************
   */
  Inlet(std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup, bool docsEnabled = false) :
    m_reader(reader),
    m_sidreRootGroup(sidreRootGroup),
    m_globalTable(std::make_shared<Table>("", "", m_reader, m_sidreRootGroup)),
    m_docWriterEnabled(docsEnabled) {}

  virtual ~Inlet() = default;

  /*!
   *****************************************************************************
   * \brief Returns the shared pointer to the Reader class.
   *
   * Provides access to the Reader class that is used to access the input deck.
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
  // Functions that define the input deck schema
  //

  /*!
   *****************************************************************************
   * \brief Add a Table to the input deck schema.
   *
   * Adds a Table to the input deck schema. Tables hold a varying amount Fields
   * defined by the user.  By default, it is not required unless marked with
   * Table::required(). This creates the Sidre Group class with the given name and
   * stores the given description.
   *
   * \param [in] name Name of the Table expected in the input deck
   * \param [in] description Description of the Table
   *
   * \return Shared pointer to the created Table
   *****************************************************************************
   */
  std::shared_ptr<Table> addTable(const std::string& name,
                                  const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Boolean Field to the input deck schema.
   *
   * Adds a Boolean Field to the input deck schema. It may or may not be required
   * to be present in the input deck. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input deck the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input deck
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addBool(const std::string& name,
                                 const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Double Field to the input deck schema.
   *
   * Adds a Double Field to the input deck schema. It may or may not be required
   * to be present in the input deck. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input deck the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input deck
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addDouble(const std::string& name,
                                   const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a Integer Field to the input deck schema.
   *
   * Adds a Integer Field to the input deck schema. It may or may not be required
   * to be present in the input deck. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input deck the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Field expected in the input deck
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Field> addInt(const std::string& name,
                                const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a String Field to the input deck schema.
   *
   * Adds a String Field to the input deck schema. It may or may not be required
   * to be present in the input deck. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input deck the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Table expected in the input deck
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
   * been actually present in the input deck and will be indicted by the return
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
   * been actually present in the input deck and will be indicted by the return
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
   * been actually present in the input deck and will be indicted by the return
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
   * been actually present in the input deck and will be indicted by the return
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
   * \brief Writes input deck documentation.
   *
   * This writes the input deck's documentation through the registered DocWriter.
   *
   *****************************************************************************
   */
  void writeDocs();

  // TODO add update value functions
private:
  axom::sidre::View* baseGet(const std::string& name);

  std::shared_ptr<Reader> m_reader;
  axom::sidre::Group* m_sidreRootGroup = nullptr;
  std::shared_ptr<Table> m_globalTable;

  std::shared_ptr<DocWriter> m_docWriter;
  bool m_docWriterEnabled = false; 
};

} // end namespace inlet
} // end namespace axom

#endif
