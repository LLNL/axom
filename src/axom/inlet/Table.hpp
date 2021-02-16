// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
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
#include "axom/inlet/Function.hpp"
#include "axom/inlet/Reader.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/VariantKey.hpp"
#include "axom/inlet/Verifiable.hpp"

#include "axom/sidre.hpp"

/*!
 *******************************************************************************
 * \brief Prototype for user-defined types wishing to define a by-value read
 * from inlet with axom::inlet::Table::get(const std::string&)
 * 
 * This is the only way of reading in a non-default-constructible type from
 * Inlet.
 * 
 * \see axom::inlet::Table::get(const std::string&)
 *******************************************************************************
 */
template <typename T>
struct FromInlet
{
  // "Poison" the base implementation so it cannot be used
  // and so it can be checked in a SFINAE context
  FromInlet() = delete;
};

namespace axom
{
namespace inlet
{
// Forward declaration for the traits

class Table;

namespace detail
{
/*!
 *******************************************************************************
 * \class is_inlet_primitive
 *
 * \brief A type trait for checking if a type is isomorphic to an Inlet primitive
 * \tparam T The type to check
 * \note An Inlet primitive is any of the following C++ types: int, double, bool,
 * std::string.
 * \see InletType
 *******************************************************************************
 */
template <typename T>
struct is_inlet_primitive
{
  using BaseType = typename std::decay<T>::type;
  static constexpr bool value = std::is_same<BaseType, bool>::value ||
    std::is_same<BaseType, int>::value || std::is_same<BaseType, double>::value ||
    std::is_same<BaseType, std::string>::value;
};

/*!
 *******************************************************************************
 * \class is_inlet_primitive
 *
 * \brief A type trait for checking if a type is isomorphic to an array of Inlet
 * primitives
 * \tparam T The type to check
 *******************************************************************************
 */
template <typename T>
struct is_inlet_primitive_array : std::false_type
{ };

// If it's an unordered map whose value type is a Inlet primitive,
// assume that it's an array
template <typename T>
struct is_inlet_primitive_array<std::unordered_map<int, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_inlet_array : std::false_type
{ };

template <typename T>
struct is_inlet_array<std::unordered_map<int, T>> : std::true_type
{ };

// Determines whether a type is a std::function
template <typename T>
struct is_std_function : std::false_type
{ };

template <typename T>
struct is_std_function<std::function<T>> : std::true_type
{ };

// Extracts the signature of a std::function
template <typename FuncType>
struct std_function_signature;

template <typename FuncType>
struct std_function_signature<std::function<FuncType>>
{
  using type = FuncType;
};

template <typename T>
struct is_inlet_dict : std::false_type
{ };

template <typename T>
struct is_inlet_dict<std::unordered_map<std::string, T>> : std::true_type
{ };

template <typename T>
struct is_inlet_dict<std::unordered_map<VariantKey, T>> : std::true_type
{ };

template <typename T>
struct is_inlet_primitive_dict : std::false_type
{ };

template <typename T>
struct is_inlet_primitive_dict<std::unordered_map<std::string, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_inlet_primitive_dict<std::unordered_map<VariantKey, T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

template <typename T>
struct is_std_vector : std::false_type
{ };

template <typename T>
struct is_std_vector<std::vector<T>> : std::true_type
{ };

template <typename T>
struct is_primitive_std_vector : std::false_type
{ };

template <typename T>
struct is_primitive_std_vector<std::vector<T>>
{
  static constexpr bool value = is_inlet_primitive<T>::value;
};

/*!
 *******************************************************************************
 * \class has_FromInlet_specialization
 *
 * \brief A type trait for checking if a type has specialized FromInlet
 * with the required T operator()(axom::inlet::Table&)
 * \tparam T The type to check
 *******************************************************************************
 */
template <typename T, typename SFINAE = void>
struct has_FromInlet_specialization : std::false_type
{ };

template <typename T>
struct has_FromInlet_specialization<
  T,
  typename std::enable_if<std::is_same<T,
                                       decltype(std::declval<FromInlet<T>&>()(
                                         std::declval<const Table&>()))>::value>::type>
  : std::true_type
{ };

/*!
 *******************************************************************************
 * \brief An overloaded utility function for converting a type to a string
 * 
 * \note Needed as std::to_string doesn't implement an identity overload
 * 
 * \param [in] idx The index to convert to string
 *******************************************************************************
 */
inline std::string indexToString(const std::string& idx) { return idx; }
/// \overload
inline std::string indexToString(const int idx) { return std::to_string(idx); }
/// \overload
inline std::string indexToString(const VariantKey& idx)
{
  return idx.type() == InletType::String ? static_cast<std::string>(idx)
                                         : indexToString(static_cast<int>(idx));
}

/*!
 *******************************************************************************
 * \brief An templated utility function for converting an index to the desired type
 * 
 * \param [in] idx The index to convert
 * \tparam From The type of the idx parameter (converting from)
 * \tparam Result The type to convert to
 *******************************************************************************
 */
template <typename Result, typename From>
inline Result toIndex(const From& idx)
{
  // By default, try a static cast
  return static_cast<Result>(idx);
}
/// \overload
template <>
inline int toIndex(const std::string& idx)
{
  int idx_as_int;
  SLIC_ERROR_IF(!checkedConvertToInt(idx, idx_as_int),
                fmt::format("[Inlet] Expected an integer, got: {0}", idx));
  return idx_as_int;
}

/*!
 *******************************************************************************
 * \brief Determines whether a variant key is convertible to another type
 * 
 * \tparam Key The type to check the validity of the conversion to
 * 
 * \note That is, returns true if the key holds an integer and Key is int, etc
 *******************************************************************************
 */
template <typename Key>
bool matchesKeyType(const VariantKey& key)
{
  // A VariantKey matches everything
  if(std::is_same<Key, VariantKey>::value)
  {
    return true;
  }
  else if(std::is_same<Key, int>::value && key.type() == InletType::Integer)
  {
    return true;
  }
  else if(std::is_same<Key, std::string>::value && key.type() == InletType::String)
  {
    return true;
  }
  return false;
}

}  // namespace detail

class Proxy;
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
class Table : public Verifiable<Table>
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
   * \param [in] reader Reference to the input file Reader class.
   * \param [in] sidreRootGroup Pointer to the already created Sidre Group.
   * \param [in] docEnabled Boolean indicating whether or not documentation
   * generation is enabled for input feck this Table instance belongs to.
   *****************************************************************************
   */
  Table(const std::string& name,
        const std::string& description,
        Reader& reader,
        axom::sidre::Group* sidreRootGroup,
        bool docEnabled = true)
    : m_name(name)
    , m_reader(reader)
    , m_sidreRootGroup(sidreRootGroup)
    , m_docEnabled(docEnabled)
  {
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

  // Tables must be move-only - delete the implicit shallow copy constructor
  Table(const Table&) = delete;
  Table(Table&&) = default;

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
  /// \overload
  const axom::sidre::Group* sidreGroup() const { return m_sidreGroup; };

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
   * \brief Add an array of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addBoolArray(const std::string& name,
                                  const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Integer Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addIntArray(const std::string& name,
                                 const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Double Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addDoubleArray(const std::string& name,
                                    const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of String Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Verifiable<Table>& addStringArray(const std::string& name,
                                    const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   *
   * \return Reference to the created array
   *****************************************************************************
   */
  Table& addStructArray(const std::string& name,
                        const std::string& description = "");
  /*!
   *****************************************************************************
   * \brief Add a dictionary of Boolean Fields to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the Field
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Verifiable<Table>& addBoolDictionary(const std::string& name,
                                       const std::string& description = "");

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
                                      const std::string& description = "");
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
                                         const std::string& description = "");

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
                                         const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a dictionary of user-defined types to the input file schema.
   *
   * \param [in] name Name of the dict
   * \param [in] description Description of the dictionary
   *
   * \return Reference to the created dictionary
   *****************************************************************************
   */
  Table& addStructDictionary(const std::string& name,
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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addBool(const std::string& name,
                            const std::string& description = "")
  {
    return addPrimitive<bool>(name, description);
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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addDouble(const std::string& name,
                              const std::string& description = "")
  {
    return addPrimitive<double>(name, description);
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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addInt(const std::string& name,
                           const std::string& description = "")
  {
    return addPrimitive<int>(name, description);
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
   * \return Reference to the created Field
   *****************************************************************************
   */
  VerifiableScalar& addString(const std::string& name,
                              const std::string& description = "")
  {
    return addPrimitive<std::string>(name, description);
  }

  /*!
   *****************************************************************************
   * \brief Get a function from the input deck
   *
   * \param [in] name         Name of the function
   * \param [in] ret_type     The return type of the function
   * \param [in] arg_types    The argument types of the function
   * \param [in] description  Description of the function
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Function
   *****************************************************************************
   */
  Verifiable<Function>& addFunction(const std::string& name,
                                    const FunctionTag ret_type,
                                    const std::vector<FunctionTag>& arg_types,
                                    const std::string& description = "",
                                    const std::string& pathOverride = "");

  /*!
   *******************************************************************************
   * \brief Returns a stored value of primitive type.
   * 
   * Retrieves a value of primitive type.
   * 
   * \param [in] name Name of the Field value to be gotten
   * \return The retrieved value
   * 
   * \tparam T The primitive type
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_primitive<T>::value, T>::type get(
    const std::string& name) const
  {
    if(!hasField(name))
    {
      const std::string msg = fmt::format(
        "[Inlet] Field with specified path"
        "does not exist: {0}",
        name);
      SLIC_ERROR(msg);
    }
    return getField(name).get<T>();
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
   * \pre Requires a specialization of FromInlet for T
   * 
   * \tparam T The user-defined type
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<
    !detail::is_inlet_primitive<T>::value && !detail::is_inlet_array<T>::value &&
      !detail::is_inlet_dict<T>::value && !detail::is_std_vector<T>::value,
    T>::type
  get(const std::string& name = "") const
  {
    static_assert(detail::has_FromInlet_specialization<T>::value,
                  "To read a user-defined type, specialize FromInlet<T>");
    FromInlet<T> from_inlet;
    if(name.empty())
    {
      return from_inlet(*this);
    }
    else
    {
      if(!hasTable(name))
      {
        std::string msg =
          fmt::format("[Inlet] Table with name '{0}' does not exist", name);
        SLIC_ERROR(msg);
      }
      return from_inlet(getTable(name));
    }
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored collection.
   * 
   * Retrieves a collection of user-defined type.
   * 
   * \return The retrieved collection
   * 
   * \tparam T The collection type, i.e., T = std::unordered_map<K, V>
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_inlet_array<T>::value ||
                            detail::is_inlet_dict<T>::value,
                          T>::type
  get() const
  {
    using Key = typename T::key_type;
    using Val = typename T::mapped_type;
    // This needs to work transparently for both references to the underlying
    // internal table and references using the same path as the data file
    if(isCollectionGroup(m_name))
    {
      return getCollection<Key, Val>();
    }
    // If the collection group is a child table, retrieve it and copy its contents
    // into the result
    else
    {
      return getTable(detail::COLLECTION_GROUP_NAME).getCollection<Key, Val>();
    }
  }

  /*!
   *******************************************************************************
   * \brief Returns a stored collection as a contiguous array.
   * 
   * Retrieves a collection of user-defined type.
   * 
   * \return The values in the retrieved collection
   * 
   * \tparam T The collection type, i.e., T = std::vector<V>
   * 
   * \note Elements in the returned array will be in ascending order by index,
   * regardless of index contiguity or base index
   *******************************************************************************
   */
  template <typename T>
  typename std::enable_if<detail::is_std_vector<T>::value, T>::type get() const
  {
    // Only allow retrieval of std::vectors from integer-keyed collections
    using Key = int;
    using Val = typename T::value_type;
    auto map = get<std::unordered_map<Key, Val>>();

    // Retrieve and sort the indices to provide consistent behavior regardless
    // of index contiguity or base index
    std::vector<Key> indices;
    indices.reserve(map.size());

    for(const auto& entry : map)
    {
      indices.push_back(entry.first);
    }
    std::sort(indices.begin(), indices.end());

    std::vector<Val> result;
    result.reserve(map.size());

    for(const Key index : indices)
    {
      // Safe to move from the map as it won't get used afterwards
      result.push_back(std::move(map[index]));
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
   *******************************************************************************
   */
  Proxy operator[](const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Set the required status of this Table.
   *
   * Set whether this Table is required, or not, to be in the input file.
   * The default behavior is to not be required.
   *
   * \param [in] isRequired Boolean value of whether Table is required
   *
   * \return Reference to this instance of Table
   *****************************************************************************
   */
  Table& required(bool isRequired = true);

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
  bool isRequired() const;

  /*!
   *****************************************************************************
   * \brief Registers the function object that will verify this Table's contents
   * during the verification stage.
   * 
   * \param [in] The function object that will be called by Table::verify().
   *****************************************************************************
  */
  Table& registerVerifier(std::function<bool(const Table&)> lambda);

  /*!
   *****************************************************************************
   * \brief This will be called by Inlet::verify to verify the contents of this
   *  Table and all child Tables/Fields of this Table.
   *****************************************************************************
  */
  bool verify() const;

  /*!
   *****************************************************************************
   * \brief Return whether a Table or Field with the given name is present in 
   * this Table's subtree.
   *
   * \return Boolean value indicating whether this Table's subtree contains a
   * Field or Table with the given name.
   *****************************************************************************
   */
  bool contains(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Returns whether this table or any of its subtables contain a non-
   * empty field
   *****************************************************************************
   */
  explicit operator bool() const;

  /*!
   *****************************************************************************
   * \return An unordered map from Field names to the child Field pointers for 
   * this Table.
   *****************************************************************************
   */
  const std::unordered_map<std::string, std::unique_ptr<Field>>& getChildFields() const;

  /*!
   *****************************************************************************
   * \return An unordered map from Table names to the child Table pointers for 
   * this Table.
   *****************************************************************************
   */
  const std::unordered_map<std::string, std::unique_ptr<Table>>& getChildTables() const;

  /*!
   *****************************************************************************
   * \return The full name of this Table.
   *****************************************************************************
   */
  std::string name() const;

  /*!
   *****************************************************************************
   * \brief Add a Field to the input file schema.
   *
   * Adds a primitive Field to the input file schema. It may or may not be required
   * to be present in the input file. This creates the Sidre Group class with the
   * given name and stores the given description. If present in the input file the
   * value is read and stored in the datastore. 
   *
   * \param [in] name Name of the Table expected in the input file
   * \param [in] description Description of the Table
   * \param [in] forArray Whether the primitive is in an array, in which
   * case the provided value should be inserted instead of the one read from
   * the input file
   * \param [in] val A provided value, will be overwritten if found at specified
   * path in input file
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Field
   *****************************************************************************
   */

  template <typename T,
            typename SFINAE =
              typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  VerifiableScalar& addPrimitive(const std::string& name,
                                 const std::string& description = "",
                                 bool forArray = false,
                                 T val = T {},
                                 const std::string& pathOverride = "");

private:
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
   * \return Reference to the created Table
   *****************************************************************************
   */
  Table& addTable(const std::string& name, const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add an array of primitive Fields to the input file schema.
   *
   * \param [in] name Name of the array
   * \param [in] description Description of the Field
   * \param [in] isDict Whether to use string-valued keys for the collection
   * \param [in] pathOverride The path within the input file to read from, if
   * different than the structure of the Sidre datastore
   *
   * \return Reference to the created Field
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE =
              typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  Verifiable<Table>& addPrimitiveArray(const std::string& name,
                                       const std::string& description = "",
                                       const bool isDict = false,
                                       const std::string& pathOverride = "");

  /*!
   *****************************************************************************
   * \brief Return whether a Table with the given name is present in this Table's subtree.
   *
   * \return Boolean value indicating whether the calling Table's subtree
   * contains a Table with the given name.
   *****************************************************************************
   */
  bool hasTable(const std::string& tableName) const;

  /*!
   *****************************************************************************
   * \brief Return whether a Field with the given name is present in this Table's
   *  subtree.
   *
   * \return Boolean value indicating whether the calling Table's subtree
   * contains a Field with the given name.
   *****************************************************************************
   */
  bool hasField(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Return whether a Function with the given name is present in this Table's
   *  subtree.
   *
   * \return Boolean value indicating whether the calling Table's subtree
   * contains a Function with the given name.
   *****************************************************************************
   */
  bool hasFunction(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Table.
   * 
   * \param [in] The string indicating the target name of the Table to be searched for.
   * 
   * \return The Table matching the target name.
   *****************************************************************************
   */
  Table& getTable(const std::string& tableName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Field.
   * 
   * \param [in] The string indicating the target name of the Field to be searched for.
   * 
   * \return The Field matching the target name.
   *****************************************************************************
   */
  Field& getField(const std::string& fieldName) const;

  /*!
   *****************************************************************************
   * \brief Retrieves the matching Function.
   * 
   * \param [in] The string indicating the target name of the Function to be searched for.
   * 
   * \return The Function matching the target name.
   *****************************************************************************
   */
  Function& getFunction(const std::string& funcName) const;

  /*!
   *****************************************************************************
   * \brief Helper method template for adding primitives
   * 
   * Adds the value at the templated type to the sidre group
   * 
   * \param [inout] sidreGroup The group to add the primitive view to
   * \param [in] lookupPath The path within the input file to read from
   * \param [in] forArray Whether the primitive is in an array, in which
   * case the provided value should be inserted instead of the one read from
   * the input file
   * \param [in] val A provided value, will be overwritten if found at specified
   * path in input file
   *
   * \return Type ID for the inserted view
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE =
              typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  axom::sidre::DataTypeId addPrimitiveHelper(axom::sidre::Group* sidreGroup,
                                             const std::string& lookupPath,
                                             bool forArray,
                                             T val);

  /*!
   *****************************************************************************
   * \brief Helper method template for adding primitives
   * 
   * Reads an array at the provided path into the provided table
   * 
   * \param [inout] table The inlet::Table to add the array to 
   * \param [in] lookupPath The path within the input file to read from
   * \param [in] isDict Whether to use string keys
   * 
   *****************************************************************************
   */
  template <typename T,
            typename SFINAE =
              typename std::enable_if<detail::is_inlet_primitive<T>::value>::type>
  void addPrimitiveArrayHelper(Table& table,
                               const std::string& lookupPath,
                               bool isDict = false);

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
   * \return The child Field matching the target name.
   *****************************************************************************
   */
  Field& addField(axom::sidre::Group* sidreGroup,
                  axom::sidre::DataTypeId type,
                  const std::string& fullName,
                  const std::string& name);

  /*!
   *****************************************************************************
   * \brief Adds the Function.
   * 
   * \param [in] The Sidre Group corresponding to the Function that will be added.
   * \param [in] func The actual callable to store
   * \param [in] The complete Table sequence for the Table this Function will be added to.
   * \param [in] The Table sequence for the Table this Function will be added to, 
   * relative to this Table.
   * 
   * \return The child Function matching the target name.
   *****************************************************************************
   */
  Function& addFunctionInternal(axom::sidre::Group* sidreGroup,
                                FunctionVariant&& func,
                                const std::string& fullName,
                                const std::string& name);

  axom::sidre::View* baseGet(const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief This is an internal helper that returns the pointer-to-member for
   * the unordered_map of children of requested type.
   * 
   * \tparam T The type of the child to search for (Field/Table/Function)
   *****************************************************************************
   */
  template <typename T>
  static std::unordered_map<std::string, std::unique_ptr<T>> Table::*getChildren();

  /*!
   *****************************************************************************
   * \brief This is an internal helper. It return whether this Table has a child 
   * with the given name and type.
   * 
   * \tparam T The type of the child to search for (Field/Table/Function)
   * \return Boolean value of whether this Table has the child.
   *****************************************************************************
   */
  template <typename T>
  bool hasChild(const std::string& childName) const;

  /*!
   *****************************************************************************
   * \brief This retrieves the child of requested name and type.
   * 
   * \param [in] The string indicating the target name of the child to be searched for.
   * 
   * \tparam T The type of the child to search for (Field/Table/Function)
   * \return The child matching the target name. If no such child is found,
   * a nullptr is returned.
   *****************************************************************************
   */
  template <typename T>
  T* getChildInternal(const std::string& childName) const;

  /*!
   *****************************************************************************
   * \brief This is an internal utility intended to be used with arrays/dicts of 
   * user-defined types that returns the indices as strings - integer indices
   * will be converted to strings
   * 
   * \param [in] trimAbsolute Whether to only return the "basename" if the path
   * is absolute, e.g., an absolute path foo/0/bar will be trimmed to "bar"
   *****************************************************************************
   */
  std::vector<VariantKey> collectionIndices(bool trimAbsolute = true) const;

  /*!
   *****************************************************************************
   * \brief This is an internal utility intended to be used with arrays of 
   * user-defined types that returns the a list of pairs, each of which contain
   * an index (a number) and a fully qualified path within the input file to
   * the array element at the corresponding index.
   * 
   * \param [in] name The name of the array object in the input file
   *****************************************************************************
   */
  std::vector<std::pair<std::string, std::string>> collectionIndicesWithPaths(
    const std::string& name) const;

  /*!
   *****************************************************************************
   * \brief Get a collection represented as an unordered map from the input file
   *****************************************************************************
   */
  template <typename Key, typename Val>
  std::unordered_map<Key, Val> getCollection() const
  {
    std::unordered_map<Key, Val> map;
    for(const auto& indexLabel : collectionIndices())
    {
      if(detail::matchesKeyType<Key>(indexLabel))
      {
        map[detail::toIndex<Key>(indexLabel)] =
          get<Val>(detail::indexToString(indexLabel));
      }
    }
    return map;
  }

  /*!
   *******************************************************************************
   * \brief Adds a group containing the indices of a collection to the calling 
   * table and a subtable for each index
   * 
   * \param [in] indices The indices to add
   * \param [in] description The optional description of the subtables
   * \tparam Key The type of the indices to add
   *******************************************************************************
   */
  template <typename Key>
  void addIndicesGroup(const std::vector<Key>& indices,
                       const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Add a collection of user-defined type to the input file schema.
   *
   * \param [in] name Name of the collection
   * \param [in] description Description of the collection
   *
   * \return Reference to the created collection
   *****************************************************************************
   */
  template <typename Key>
  Table& addStructCollection(const std::string& name,
                             const std::string& description = "");

  /*!
   *****************************************************************************
   * \brief Returns true if the calling object is part of a struct collection,
   * i.e., an array or dictionary of user-defined type
   *****************************************************************************
   */
  bool isStructCollection() const
  {
    return m_sidreGroup->hasView(detail::STRUCT_COLLECTION_FLAG);
  }

  /*!
   *****************************************************************************
   * \brief Calls a function on the subtables corresponding to the elements
   * of the collection held by this table
   * 
   * \param [in] func The function to apply to individual collection elements
   * 
   * \pre The calling table must be a struct collection, i.e., isStructCollection()
   * returns true
   * 
   * \pre The function must accept a single argument of type Table&
   * 
   *****************************************************************************
   */
  template <typename Func>
  void forEachCollectionElement(Func&& func) const;

  std::string m_name;
  Reader& m_reader;
  // Inlet's Root Sidre Group
  axom::sidre::Group* m_sidreRootGroup;
  // This Table's Sidre Group
  axom::sidre::Group* m_sidreGroup;
  bool m_docEnabled;
  std::unordered_map<std::string, std::unique_ptr<Table>> m_tableChildren;
  std::unordered_map<std::string, std::unique_ptr<Field>> m_fieldChildren;
  std::unordered_map<std::string, std::unique_ptr<Function>> m_functionChildren;
  std::function<bool(const Table&)> m_verifier;

  // Used for ownership only - need to take ownership of these so children
  // and their aggregates have identical lifetime
  std::vector<AggregateVerifiable<Table>> m_aggregate_tables;
  std::vector<AggregateField> m_aggregate_fields;
  std::vector<AggregateVerifiable<Function>> m_aggregate_funcs;

  // Used when the calling Table is a struct collection within a struct collection
  // Need to delegate schema-defining calls (add*) to the elements of the nested
  // collection
  std::vector<std::reference_wrapper<Table>> m_nested_aggregates;
};

}  // namespace inlet
}  // namespace axom

#endif
