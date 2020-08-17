// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Table.hpp
 *
 * \brief This file contains the class definition of Inlet's Table class.
 *******************************************************************************
 */

#ifndef INLET_TABLE_HPP
#define INLET_TABLE_HPP

#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include <tuple>

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/SchemaCreator.hpp"

#include "axom/sidre.hpp"

namespace axom
{
namespace inlet
{

/*!
 *******************************************************************************
 * \class Table
 *
 * \brief Provides functions to help define how individual Table and Field
 *        variables in an input deck are expected to behave.  It also holds the
 *        Sidre Group to the individual Table.
 *
 * \see Inlet Field
 *******************************************************************************
 */
class Table : public SchemaCreator, public std::enable_shared_from_this<Table>
{
public:
  /*!
   *****************************************************************************
   * \brief Constructor for the Table class.
   *
   * This class provides functions to define the behavior of the Table
   * data already read and stored in the given Sidre Group. This creates
   * the necessary Sidre Group's and views to store the given name and description.
   *
   * \param [in] name Name of the Table expected in the input deck
   * \param [in] description Description of the Table
   * \param [in] reader Shared pointer to the input deck Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for Input Deck this Table instance belongs to.
   *****************************************************************************
   */
  Table(const std::string& name,
        const std::string& description,
        std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true) : 
    m_name(name),
    m_reader(reader),
    m_sidreRootGroup(sidreRootGroup),
    m_docEnabled(docEnabled)
    {
      SLIC_ASSERT_MSG(m_reader, "Inlet's Reader class not valid");
      SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr, "Inlet's Sidre Datastore class not set");

      if (m_name == "")
      {
        m_sidreGroup = m_sidreRootGroup;
      }
      else
      {
        if (!m_sidreRootGroup->hasGroup(name))
        {
          m_sidreGroup = m_sidreRootGroup->createGroup(name);
        }
        else
        {
          m_sidreGroup = m_sidreRootGroup->getGroup(name);
        }
      }

      if(description != "")
      {
        if (m_sidreGroup->hasView("description"))
        {
          //TODO: warn user?
          m_sidreGroup->destroyViewAndData("description");
        }
        m_sidreGroup->createViewString("description", description);
      }
    }

  virtual ~Table() = default;

  /*!
   *****************************************************************************
   * \brief Returns pointer to the Sidre Group class for this Table.
   *
   * Provides access to the Sidre Group class that holds all the stored
   * information for this Table class.
   *
   * \return Pointer to the Sidre Group for this Table
   *****************************************************************************
   */
  axom::sidre::Group* sidreGroup() { return m_sidreGroup; };

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
                                  const std::string& description="");

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
                                 const std::string& description="");

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
                                   const std::string& description="");

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
                                const std::string& description="");

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
                                   const std::string& description="");

  /*!
   *****************************************************************************
   * \brief Set the required status of this Table.
   *
   * Set whether this Table is required, or not, to be in the input deck.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Table is required
   *
   * \return Shared pointer to this instance of Table
   *****************************************************************************
   */
  std::shared_ptr<Table> required(bool isRequired);

  /*!
   *****************************************************************************
   * \brief Return the required status of this Table.
   *
   * Return that this Table is required, or not, to be in the input deck.
   * The default behavior is to not be required.
   *
   * \return Boolean value of whether this Table is required
   *****************************************************************************
   */
  bool required();

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this Table's contents
   * during the verification stage.
   * 
   * \param [in] The function object that will be called by Table::verify().
   *****************************************************************************
  */
  std::shared_ptr<Table> registerVerifier(std::function<bool()> lambda);

   /*!
   *****************************************************************************
   * \brief This will be called by Inlet::verify to verify the contents of this
   *  Table and all child Tables/Fields of this Table.
   *****************************************************************************
  */
  bool verify();

  /*!
   *****************************************************************************
   * \brief Return whether this Table has a child Field with the given name.
   *
   * \return Boolean value of whether this Table has the child Field.
   *****************************************************************************
   */
  bool hasField(const std::string& fieldName);

  /*!
   *****************************************************************************
   * \brief Return whether this Table has a child Table with the given name.
   *
   * \return Boolean value of whether this Table has the child Table.
   *****************************************************************************
   */
  bool hasTable(const std::string& tableName);

  /*!
   *****************************************************************************
   * \return An unordered map from Field names to the child Field pointers for 
   * this Table.
   *****************************************************************************
   */
  std::unordered_map<std::string, std::shared_ptr<Field>> getChildFields();

  /*!
   *****************************************************************************
   * \return An unordered map from Table names to the child Table pointers for 
   * this Table.
   *****************************************************************************
   */
  std::unordered_map<std::string, std::shared_ptr<Table>> getChildTables();

  /*!
   *****************************************************************************
   * \return The full name of this Table.
   *****************************************************************************
   */
  std::string name();

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Table.
   * 
   * \param [in] The string indicating the target name of the Table to be searched for.
   * 
   * \return The child Table matching the target name. If no such Table is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Table> getTable(const std::string& tableName);

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
  std::shared_ptr<Field> getField(const std::string& fieldName);

private:
  /*!
   *****************************************************************************
   * \brief Creates the basic Sidre Group for this Table and stores the given
   *        information
   *
   * \return Pointer to the created Sidre Group for this Table
   *****************************************************************************
   */
  axom::sidre::Group* createSidreGroup(const std::string& name,
                                   const std::string& description);

  /*!
   *****************************************************************************
   * \brief Adds the Field.
   * 
   * \param [in] The Sidre Group corresponding to the Field that will be added.
   * \param [in] The type ID
   * \param [in] The complete Table sequence for the Table this Field will be added to.
   * \param [in] The Table sequence for the Table this Field will be added to, 
   * relative to this Table.
   * 
   * \return The child Field matching the target name. If no such Field is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Field> addField(axom::sidre::Group* sidreGroup, 
                                  axom::sidre::DataTypeId type, 
                                  const std::string& fullName,
                                  const std::string& name);

  std::string m_name;
  std::shared_ptr<Reader> m_reader;
  // Inlet's Root Sidre Group
  axom::sidre::Group* m_sidreRootGroup;
  // This Table's Sidre Group
  axom::sidre::Group* m_sidreGroup;
  bool m_docEnabled;
  std::unordered_map<std::string, std::shared_ptr<Table>> m_tableChildren;
  std::unordered_map<std::string, std::shared_ptr<Field>> m_fieldChildren;
  std::function<bool()> m_verifier;
};

} // end namespace inlet
} // end namespace axom

#endif
