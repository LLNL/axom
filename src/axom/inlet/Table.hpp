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
#include <type_traits>

#include "fmt/fmt.hpp"

#include "axom/inlet/Field.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/SchemaCreator.hpp"

#include "axom/sidre.hpp"

/*!
 *******************************************************************************
 * \brief Prototype for user-defined types wishing to define a by-value read
 * from inlet with axom::inlet::Table::get(const std::string&)
 * 
 * This is the only way of reading in a non-default-constructible type from
 * Inlet.  Usage of axom::inlet::Table::get(const std::string&) without
 * specializing this function will result in a linker error.
 * 
 * \see axom::inlet::Table::get(const std::string&)
 * 
 * \param [in] base The Inlet table representing the root of the object
 * in which all object fields are children
 *******************************************************************************
 */
template <typename T>
T from_inlet(axom::inlet::Table& base);

namespace axom
{
namespace inlet
{
template <typename T>
struct is_lua_primitive
{
  using BaseType = typename std::decay<T>::type;
  static constexpr bool value = std::is_same<BaseType, bool>::value ||
    std::is_same<BaseType, int>::value || std::is_same<BaseType, double>::value ||
    std::is_same<BaseType, std::string>::value;
};

// By default it's false
template <typename T>
struct is_lua_primitive_array
{
  static constexpr bool value = false;
};

// If it's an unordered map whose value type is a lua primitive,
// assume that it's an array
template <typename T>
struct is_lua_primitive_array<std::unordered_map<int, T>>
{
  static constexpr bool value = is_lua_primitive<T>::value;
};

/*!
 *******************************************************************************
 * \class Proxy
 *
 * \brief Provides a uniform interface for access and conversion to primitive
 * and user-defined types
 *
 * \see Inlet Field
 * \see Inlet Table
 *******************************************************************************
 */
class Proxy
{
public:
  Proxy() = default;
  Proxy(Table& table) : m_table(&table) { }
  Proxy(Field& field) : m_field(&field) { }

  template <typename T>
  operator T()
  {
    return get<T>();
  }

  // user-defined types
  template <typename T>
  typename std::enable_if<!is_lua_primitive<T>::value, T>::type get();

  // primitives
  template <typename T>
  typename std::enable_if<is_lua_primitive<T>::value, T>::type get();

private:
  Table* m_table = nullptr;
  Field* m_field = nullptr;
};

/*!
 *******************************************************************************
 * \class Table
 *
 * \brief Provides functions to help define how individual Table and Field
 *        variables in an input file are expected to behave.  It also holds the
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
   * \param [in] name Name of the Table expected in the input file
   * \param [in] description Description of the Table
   * \param [in] reader Shared pointer to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for input feck this Table instance belongs to.
   *****************************************************************************
   */
  Table(const std::string& name,
        const std::string& description,
        std::shared_ptr<Reader> reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true)
    : m_name(name)
    , m_reader(reader)
    , m_sidreRootGroup(sidreRootGroup)
    , m_docEnabled(docEnabled)
  {
    SLIC_ASSERT_MSG(m_reader, "Inlet's Reader class not valid");
    SLIC_ASSERT_MSG(m_sidreRootGroup != nullptr,
                    "Inlet's Sidre Datastore class not set");

    if(m_name == "")
    {
      m_sidreGroup = m_sidreRootGroup;
    }
    else
    {
      if(!m_sidreRootGroup->hasGroup(name))
      {
        m_sidreGroup = m_sidreRootGroup->createGroup(name);
        m_sidreGroup->createViewString("InletType", "Table");
      }
      else
      {
        m_sidreGroup = m_sidreRootGroup->getGroup(name);
      }
    }

    if(description != "")
    {
      if(m_sidreGroup->hasView("description"))
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
   * \brief Add an array of Boolean Fields to the input deck schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Shared pointer to the created Field
   *****************************************************************************
   */
  std::shared_ptr<Table> addBoolArray(const std::string& name,
                                      const std::string& description = "");

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
                                     const std::string& description = "");

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
                                        const std::string& description = "");

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
                                        const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Get a boolean array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, bool>& map);

  /*!
   *****************************************************************************
   * \brief Get a int array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, int>& map);

  /*!
   *****************************************************************************
   * \brief Get a double array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, double>& map);

  /*!
   *****************************************************************************
   * \brief Get a string array represented as an unordered map from the input deck
   *
   * \param [out] map Unordered map to be populated with array contents
   *
   * \return Whether or not the array was found
   *****************************************************************************
   */
  bool getArray(std::unordered_map<int, std::string>& map);

  /*!
   *****************************************************************************
   * \brief Add a Boolean Field to the input deck schema.
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
                                 const std::string& description = "")
  {
    return addBoolHelper(name, description);
  }

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
                                   const std::string& description = "")
  {
    return addDoubleHelper(name, description);
  }

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
                                const std::string& description = "")
  {
    return addIntHelper(name, description);
  }
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
                                   const std::string& description = "")
  {
    return addStringHelper(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Gets a value of primitive type out of the table with an out-param.
   *
   * Retrieves the Field value out of the table.  This Field may not have
   * been actually present in the input file and will be indicated by the return
   * value. 
   *
   * \param [in] name Name of the Field value to be gotten
   * \param [out] value Value to be filled
   *
   * \return True if the value was found in the table
   * \tparam T The type of the value to retrieve
   *****************************************************************************
   */
  template <typename T>
  typename std::enable_if<is_lua_primitive<T>::value, bool>::type get(
    const std::string& name,
    T& value)
  {
    bool found = false;
    if(hasField(name))
    {
      found = getField(name)->get(value);
    }
    return found;
  }

  /*!
 *******************************************************************************
 * \brief Gets a value of user-defined type out of the table with an out-param.
 * 
 * Retrieves a value of user-defined type.
 * 
 * \param [in] name The name of the subtable representing the root of the object
 * If the empty string is passed, the calling table is used as the root
 * \param [out] value Value to be filled
 * \return True if the value was found in the table
 * \tparam T The user-defined type to retrieve
 * \pre If T is a user-defined type, requires a function
 * \code{.cpp}
 * void from_inlet(axom::inlet::Table&, T&);
 * \endcode
 * to be defined
 *******************************************************************************
 */
  template <typename T>
  typename std::enable_if<!is_lua_primitive<T>::value, bool>::type get(
    const std::string& name,
    T& value)
  {
    bool found = false;
    if(name.empty())
    {
      found = true;
      from_inlet(*this, value);
    }
    else if(hasTable(name))
    {
      found = true;
      from_inlet(*getTable(name), value);
    }
    return found;
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored value of primitive type.
   * 
   * Retrieves a value of primitive type.
   * 
   * \param [in] name Name of the Field value to be gotten
   * \return The retrieved value
   * \exception std::out_of_range If the requested field does not exist or does not 
   * contain the requested type
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<is_lua_primitive<T>::value, T>::type get(
    const std::string& name)
  {
    if(!hasField(name))
    {
      const std::string msg = fmt::format(
        "[Inlet] Field with specified path"
        "does not exist: {0}",
        name);
      throw std::out_of_range(msg);
    }
    return getField(name)->get<T>();
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored value of user-defined type.
   * 
   * Retrieves a value of user-defined type.
   * 
   * \param [in] name The name of the subtable representing the root of the object
   * If nothing is passed, the calling table is interpreted as the roof of the object
   * \return The retrieved value
   * \pre Requires a specialization of T from_inlet<T>(axom::inlet::Table&)
   * \exception std::out_of_range If the requested subtable does not exist
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<!is_lua_primitive<T>::value &&
                            !is_lua_primitive_array<T>::value,
                          T>::type
  get(const std::string& name = "")
  {
    if(name.empty())
    {
      return from_inlet<T>(*this);
    }
    else
    {
      if(!hasTable(name))
      {
        std::string msg =
          fmt::format("[Inlet] Table with name {0} does not exist", name);
        throw std::out_of_range(msg);
      }
      return from_inlet<T>(*getTable(name));
    }
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored array of primitive types.
   * 
   * Retrieves a value of user-defined type.
   * 
   * \return The retrieved array
   * \exception std::out_of_range If the calling table does not contain an array
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<is_lua_primitive_array<T>::value, T>::type get()
  {
    T result;
    if(!getTable("_inlet_array")->getArray(result))
    {
      throw std::out_of_range(
        "[Inlet] Table does not contain a valid array of requested type");
    }
    return result;
  }

  /*!
   *******************************************************************************
   * \brief Obtains a proxy view into the table for either a Field/Table subobject
   * 
   * Returns a reference via a lightweight proxy object to the element in the 
   * datastore at the index specified by the name.  This can be a field 
   * or a table.
   * 
   * \param [in] name The name of the subobject
   * \return The retrieved array
   * \exception std::out_of_range If the calling table contains both a field and 
   * a subtable with the specified name, or if the calling table does not have 
   * either a field or a subtable with the specified name
   *******************************************************************************
   */
  Proxy operator[](const std::string& name);

  /*!
   *****************************************************************************
   * \brief Set the required status of this Table.
   *
   * Set whether this Table is required, or not, to be in the input file.
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
   * Return that this Table is required, or not, to be in the input file.
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
   * \brief Return whether a Table with the given name is present in this Table's subtree.
   *
   * \return Boolean value indicating whether this Table's subtree contains this Table.
   *****************************************************************************
   */
  bool hasTable(const std::string& tableName);

  /*!
   *****************************************************************************
   * \brief Return whether a Field with the given name is present in this Table's
   *  subtree.
   *
   * \return Boolean value indicating whether this Table's subtree contains this Field.
   *****************************************************************************
   */
  bool hasField(const std::string& fieldName);

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
   * \return The Table matching the target name. If no such Table is found,
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
   * \return The Field matching the target name. If no such Field is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Field> getField(const std::string& fieldName);

private:
  std::shared_ptr<Field> addBoolHelper(const std::string& name,
                                       const std::string& description = "",
                                       bool forArray = false,
                                       bool num = 0);
  std::shared_ptr<Field> addIntHelper(const std::string& name,
                                      const std::string& description = "",
                                      bool forArray = false,
                                      int num = 0);
  std::shared_ptr<Field> addDoubleHelper(const std::string& name,
                                         const std::string& description = "",
                                         bool forArray = false,
                                         double num = 0);
  std::shared_ptr<Field> addStringHelper(const std::string& name,
                                         const std::string& description = "",
                                         bool forArray = false,
                                         const std::string& str = "");
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

  /*!
   *****************************************************************************
   * \brief This is the internal implementation of getTable. It retrieves the matching Table.
   * 
   * \param [in] The string indicating the target name of the Table to be searched for.
   * 
   * \return The Table matching the target name. If no such Table is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Table> getTableInternal(const std::string& tableName);

  /*!
   *****************************************************************************
   * \brief This is the internal implementation of getField. It retrieves the matching Field.
   * 
   * \param [in] The string indicating the target name of the Field to be searched for.
   * 
   * \return The Field matching the target name. If no such Table is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  std::shared_ptr<Field> getFieldInternal(const std::string& fieldName);

  /*!
   *****************************************************************************
   * \brief This is an internal helper. It returns whether this Table has a child 
   * Table with the given name.
   *
   * \return Boolean value of whether this Table has the child Table.
   *****************************************************************************
   */
  bool hasChildTable(const std::string& tableName);

  /*!
   *****************************************************************************
   * \brief This is an internal helper. It return whether this Table has a child 
   * Field with the given name.
   *
   * \return Boolean value of whether this Table has the child Field.
   *****************************************************************************
   */
  bool hasChildField(const std::string& fieldName);

  axom::sidre::View* baseGet(const std::string& name);

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

template <typename T>
typename std::enable_if<is_lua_primitive<T>::value, T>::type Proxy::get()
{
  SLIC_ASSERT_MSG(
    m_field != nullptr,
    "[Inlet] Tried to read a primitive type from a Proxy containing a table");
  T result;
  bool found = m_field->get(result);
  SLIC_ASSERT_MSG(
    found,
    "[Inlet] Failed to read a primitive type from a Field-containing Proxy");
  return result;
}

template <typename T>
typename std::enable_if<!is_lua_primitive<T>::value, T>::type Proxy::get()
{
  SLIC_ASSERT_MSG(m_table != nullptr,
                  "[Inlet] Tried to read a user-defined type from a Proxy "
                  "containing a single field");
  return m_table->get<T>();
}

}  // end namespace inlet
}  // end namespace axom

#endif
