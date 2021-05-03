// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MINT_ConnectivityArray_indirection_HPP_
#define MINT_ConnectivityArray_indirection_HPP_

#include <iostream>

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"
#include "axom/core/Array.hpp"
#include "axom/mint/mesh/CellTypes.hpp"
#include "axom/mint/config.hpp"
#include "axom/mint/mesh/internal/ConnectivityArrayHelpers.hpp"
#include "axom/slic/interface/slic.hpp"

#ifdef AXOM_MINT_USE_SIDRE
  #include "axom/sidre/core/sidre.hpp"
#endif

#include <cstring>

namespace axom
{
namespace mint
{
/*!
 * \class ConnectivityArray
 *
 * \brief Provides an interface for general mesh connectivity.
 *
 *  In this specialized ConnectivityArray it is assumed that each ID is of
 *  the same type but has a different number of values, as such an offset array
 *  is used to store the starting location of the values corresponding to each
 *  ID.
 *
 * \see ConnectivityArray_indirection.hpp
 * \see ConnectivityArray_typed_indirection.hpp
 * \see ConnectivityArray_internal.hpp
 */
template <>
class ConnectivityArray<INDIRECTION>
{
public:
  /*!
   * \brief Default constructor. Disabled.
   */
  ConnectivityArray() = delete;

  /// \name Native Storage ConnectivityArray Constructors
  /// @{

  /*!
   * \brief Constructs an empty ConnectivityArray instance.
   *
   * \param [in] cell_type the cell type associated with the IDs.
   * \param [in] ID_capacity the number of IDs to allocate space for.
   * \param [in] value_capacity the number of values to allocate space for.
   *
   * \pre cell_type != UNDEFINED_CELL && cell_type < NUM_CELL_TYPES
   *
   * \post getIDCapacity() >= getNumberOfIDs()
   * \post getValueCapacity() >= getNumberOfValues()
   * \post getNumberOfIDs() == 0
   * \post getNumberOfValues() == 0
   * \post getIDType() == cell_type
   */
  ConnectivityArray(CellType cell_type,
                    IndexType ID_capacity = USE_DEFAULT,
                    IndexType value_capacity = USE_DEFAULT)
    : m_cell_type(cell_type)
    , m_values(nullptr)
    , m_offsets(new Array<IndexType>(
        axom::internal::ZERO,
        1,
        (ID_capacity == USE_DEFAULT) ? USE_DEFAULT : ID_capacity + 1))
  {
    IndexType new_value_capacity =
      internal::calcValueCapacity(0, getIDCapacity(), 0, value_capacity);
    m_values = new Array<IndexType>(axom::internal::ZERO, 1, new_value_capacity);

    m_offsets->append(0);
  }

  /// @}

  /// \name External Storage ConnectivityArray Constructors
  /// @{

  /*!
   * \brief External constructor which creates a ConnectivityArray instance to
   *  wrap the given pointers.
   *
   * \param [in] cell_type the cell type associated with the IDs.
   * \param [in] n_IDs the number of IDs.
   * \param [in] values the array of values, of length at least value_capacity.
   * \param [in] offsets the offsets array, of length at least ID_capacity + 1.
   * \param [in] ID_capacity the number of IDs able to be stored in the
   *  offsets array. If not specified the capacity is set to n_IDs.
   * \param [in] value_capacity the number of values able to be stored in the
   *  values array. if not specified the capacity is set to the number of values
   *  which is calculated using the offsets array.
   *
   * \note the offset of ID 0 must be 0.
   * \note the total number of values is stored in offsets[ n_IDs ]
   *
   * \pre cell_type != UNDEFINED_CELL && cell_type < NUM_CELL_TYPES
   * \pre n_IDs >= 0
   * \pre values != nullptr
   * \pre offsets != nullptr
   *
   * \post getIDCapacity() >= getNumberOfIDs()
   * \post getValueCapacity() >= getNumberOfValues()
   * \post getNumberOfIDs() == n_IDs
   * \post getNumberOfValues() == offsets[ n_IDs ]
   * \post getIDType() == cell_type
   */
  ConnectivityArray(CellType cell_type,
                    IndexType n_IDs,
                    IndexType* values,
                    IndexType* offsets,
                    IndexType ID_capacity = USE_DEFAULT,
                    IndexType value_capacity = USE_DEFAULT)
    : m_cell_type(cell_type)
    , m_values(nullptr)
    , m_offsets(nullptr)
  {
    SLIC_ERROR_IF(n_IDs < 0,
                  "Number of IDs must be positive, not " << n_IDs << ".");
    m_offsets = new Array<IndexType>(
      offsets,
      n_IDs + 1,
      1,
      (ID_capacity == USE_DEFAULT) ? USE_DEFAULT : ID_capacity + 1);

    if(n_IDs == 0)
    {
      (*m_offsets)[0] = 0;
    }
    SLIC_ERROR_IF((*m_offsets)[0] != 0,
                  "Invalid offsets array. "
                    << "Expected item 0 to be 0 not " << (*m_offsets)[0] << ".");

    IndexType n_values = (*m_offsets)[n_IDs];
    m_values = new Array<IndexType>(values, n_values, 1, value_capacity);
  }

  /// @}

  /// \name Sidre Storage ConnectivityArray Constructors
  /// @{

#ifdef AXOM_MINT_USE_SIDRE

  /*!
   * \brief Creates a ConnectivityArray instance from a sidre::Group which
   *  already has data.
   *
   * \param [in] group the sidre::Group to create the ConnectivityArray from.
   *
   * \note the given Group must conform to a single Blueprint topology.
   *
   * \pre group != nullptr
   *
   * \post getIDCapacity() >= getNumberOfIDs()
   * \post getValueCapacity() >= getNumberOfValues()
   */
  ConnectivityArray(sidre::Group* group)
    : m_cell_type(UNDEFINED_CELL)
    , m_values(nullptr)
    , m_offsets(nullptr)
  {
    m_cell_type = internal::initializeFromGroup(group, &m_values, &m_offsets);

    SLIC_ERROR_IF(m_values->numComponents() != 1,
                  "values array must have only 1 component not "
                    << m_values->numComponents() << ".");
  }

  /*!
   * \brief Creates an empty ConnectivityArray instance from an empty
   *  sidre::Group.
   *
   * \param [in] cell_type the cell type associated with the IDs.
   * \param [in] group the sidre::Group to create the ConnectivityArray from.
   * \param [in] coordset the name of the Blueprint coordinate set to associate
   *  this ConnectivityArray with.
   * \param [in] ID_capacity the number of IDs to allocate space for.
   * \param [in] value_capacity the number of values to allocate space for.
   *
   * \pre cell_type != UNDEFINED_CELL && cell_type < NUM_CELL_TYPES
   * \pre group != nullptr
   * \pre group->getNumGroups() == group->getNumViews() == 0
   *
   * \post getIDCapacity() >= getNumberOfIDs()
   * \post getValueCapacity() >= getNumberOfValues()
   * \post getIDType() == cell_type
   */
  ConnectivityArray(CellType cell_type,
                    sidre::Group* group,
                    const std::string& coordset,
                    IndexType ID_capacity = USE_DEFAULT,
                    IndexType value_capacity = USE_DEFAULT)
    : m_cell_type(cell_type)
    , m_values(nullptr)
    , m_offsets(nullptr)
  {
    bool create_offsets = true;
    internal::initializeGroup(group, coordset, cell_type, create_offsets);

    sidre::Group* elems_group = group->getGroup("elements");
    SLIC_ASSERT(elems_group != nullptr);

    sidre::View* offsets_view = elems_group->getView("offsets");
    m_offsets = new sidre::Array<IndexType>(
      offsets_view,
      0,
      1,
      (ID_capacity == USE_DEFAULT) ? USE_DEFAULT : ID_capacity + 1);
    SLIC_ASSERT(m_offsets != nullptr);
    m_offsets->append(0);

    IndexType new_value_capacity =
      internal::calcValueCapacity(0, getIDCapacity(), 0, value_capacity);
    sidre::View* connec_view = elems_group->getView("connectivity");
    m_values = new sidre::Array<IndexType>(connec_view, 0, 1, new_value_capacity);
    SLIC_ASSERT(m_values != nullptr);
  }

#endif

  /// @}

  /*!
   * \brief Destructor, free's the allocated arrays.
   */
  ~ConnectivityArray()
  {
    if(m_values != nullptr)
    {
      delete m_values;
    }
    if(m_offsets != nullptr)
    {
      delete m_offsets;
    }
    m_values = nullptr;
    m_offsets = nullptr;
  }

  /// \name Attribute get/set Methods
  /// @{

  /*!
   * \brief Returns the total number of IDs.
   */
  IndexType getNumberOfIDs() const { return m_offsets->size() - 1; }

  /*!
   * \brief Returns the number of IDs available for storage without resizing.
   */
  IndexType getIDCapacity() const { return m_offsets->capacity() - 1; }

  /*!
   * \brief Returns the number of values in this ConnectivityArray instance.
   */
  IndexType getNumberOfValues() const { return m_values->size(); }

  /*!
   * \brief Returns the number of values available for storage without resizing.
   */
  IndexType getValueCapacity() const { return m_values->capacity(); }

  /*!
   * \brief Resize the space to hold the specified number of IDs.
   *
   * \param [in] ID_size the number of IDs to resize the space for.
   * \param [in] value_size the number of values per ID to resize the space for.
   *
   * \note if value_size is not specified, then, if this ConnectivityArray is
   *  empty, space is allocated for MAX_CELL_NODES values for each ID. Otherwise,
   *  space is allocated based on the average number of values per ID.
   */
  void resize(IndexType ID_size, IndexType value_size = USE_DEFAULT)
  {
    m_offsets->resize(ID_size + 1);
    IndexType newValueSize = internal::calcValueCapacity(getNumberOfIDs(),
                                                         getIDCapacity(),
                                                         getNumberOfValues(),
                                                         value_size);
    m_values->resize(newValueSize);
  }

  /*!
   * \brief Reserve space for IDs and values.
   *
   * \param [in] ID_capacity the number of IDs to reserve space for.
   * \param [in] value_capacity the number of values to reserve space for.
   *
   * \note if value_capacity is not specified then if this ConnectivityArray is
   *  empty then MAX_CELL_NODES values are reserved for each ID. Otherwise the
   *  average number of values per ID are reserved for each ID.
   *
   * \post getIDCapacity() >= ID_capacity
   * \post getValueCapacity() >= value_capacity
   */
  void reserve(IndexType ID_capacity, IndexType value_capacity = USE_DEFAULT)
  {
    SLIC_ERROR_IF(isExternal() && ID_capacity > m_values->capacity(),
                  "cannot exceed initial capacity of external buffer!");

    m_offsets->reserve(ID_capacity + 1);
    IndexType new_value_capacity =
      internal::calcValueCapacity(getNumberOfIDs(),
                                  getIDCapacity(),
                                  getNumberOfValues(),
                                  value_capacity);
    m_values->reserve(new_value_capacity);
  }

  /*!
   * \brief Shrink the offsets and values arrays so that there is no extra
   *  capacity.
   *
   * \post getIDCapacity() == getNumberOfIDs()
   * \post getValueCapacity() == getNumberOfValues()
   */
  void shrink()
  {
    m_values->shrink();
    m_offsets->shrink();
  }

  /*!
   * \brief Get the resize ratio.
   */
  double getResizeRatio() const { return m_values->getResizeRatio(); }

  /*!
   * \brief Set the resize ratio.
   *
   * \param [in] ratio the new resize ratio.
   *
   * \post getResizeRatio() == ratio
   */
  void setResizeRatio(double ratio)
  {
    m_values->setResizeRatio(ratio);
    m_offsets->setResizeRatio(ratio);
  }

  /*!
   * \brief Checks if this CellConnecitivity instance has a variable number of
   *  values per ID.
   * \return true.
   */
  bool hasVariableValuesPerID() const { return true; }

  /*!
   * \brief Return true if this ConnectivityArray instance is empty.
   */
  bool empty() const { return m_values->empty(); }

  /*!
   * \brief Return true iff constructed via the external constructor.
   */
  bool isExternal() const
  {
    bool consistent = true;
    bool is_external = m_values->isExternal();
    consistent &= is_external == m_offsets->isExternal();

    SLIC_WARNING_IF(!consistent, "External state not consistent.");
    return is_external;
  }

  /*!
   * \brief Return true iff constructed via the sidre constructors.
   */
  bool isInSidre() const
  {
    bool consistent = true;
    bool is_in_sidre = m_values->isInSidre();
    consistent &= is_in_sidre == m_offsets->isInSidre();

    SLIC_WARNING_IF(!consistent, "External state not consistent.");
    return is_in_sidre;
  }

  /*
   * \brief Return a const pointer to the sidre::Group that holds the data
   *  or nullptr if the data is not in sidre.
   */
#ifdef AXOM_MINT_USE_SIDRE
  const sidre::Group* getGroup() const
  {
    if(!isInSidre())
    {
      return nullptr;
    }

    return static_cast<sidre::Array<IndexType>*>(m_values)
      ->getView()
      ->getOwningGroup()
      ->getParent();
  }
#endif

  /// @}

  /// \name Data Access Methods
  /// @{

  /*!
   * \brief Returns the number of values for the given ID.
   *
   * \param [in] ID the ID in question.
   */
  IndexType getNumberOfValuesForID(IndexType ID) const
  {
    SLIC_ASSERT((ID >= 0) && (ID < getNumberOfIDs()));
    return (*m_offsets)[ID + 1] - (*m_offsets)[ID];
  }

  /*!
   * \brief Returns the cell type of the given ID.
   *
   * \param [in] ID not used, does not need to be specified.
   */
  CellType getIDType(IndexType AXOM_NOT_USED(id) = 0) const
  {
    return m_cell_type;
  }

  /*!
   * \brief Access operator for the values of the given ID.
   *
   * \param [in] ID the ID in question.
   *
   * \return pointer to the values of the given ID, of length at least
   *  getNumberOfValuesForID( ID ).
   *
   * \pre ID >= 0 && ID < getNumberOfIDs()
   * \post cell_ptr != nullptr.
   */
  /// @{

  IndexType* operator[](IndexType ID)
  {
    SLIC_ASSERT((ID >= 0) && (ID < getNumberOfIDs()));
    return m_values->getData() + (*m_offsets)[ID];
  }

  const IndexType* operator[](IndexType ID) const
  {
    SLIC_ASSERT((ID >= 0) && (ID < getNumberOfIDs()));
    return m_values->getData() + (*m_offsets)[ID];
  }

  /// @}

  /*!
   * \brief Returns a pointer to the values array, of length
   *  getNumberOfValues().
   */
  /// @{

  IndexType* getValuePtr() { return m_values->getData(); }

  const IndexType* getValuePtr() const { return m_values->getData(); }

  /// @}

  /*!
   * \brief Returns a pointer to the offsets array, of length
   *  getNumberOfIDs() + 1.
   */
  /// @{

  IndexType* getOffsetPtr() { return m_offsets->getData(); }

  const IndexType* getOffsetPtr() const { return m_offsets->getData(); }

  /// @}

  /*!
   * \brief Returns a pointer to the types array. Since this version of the
   *  ConnectivityArray does not have a types array this function returns a
   *  null pointer.
   */
  /// @{

  CellType* getTypePtr() { return nullptr; }

  const CellType* getTypePtr() const { return nullptr; }

  /// @}

  /*!
   * \brief Append a ID.
   *
   * \param [in] values pointer to the values to append, of length at least
   *  n_values.
   * \param [in] n_values the number of values corresponding to the new ID.
   * \param [in] type not used, does not need to be specified.
   *
   * \pre values != nullptr
   */
  void append(const IndexType* values,
              IndexType n_values,
              CellType AXOM_NOT_USED(type) = UNDEFINED_CELL)
  {
    SLIC_ASSERT(values != nullptr);
    m_values->append(values, n_values);
    m_offsets->append(getNumberOfValues());
  }

  /*!
   * \brief Append multiple IDs.
   *
   * \param [in] values pointer to the values to append.
   * \param [in] n_IDs the number of IDs to append.
   * \param [in] the offsets array of length n_IDs + 1.
   * \param [in] types not used, does not need to be specified.
   *
   * \note The number of values to append is given by
   *  offsets[n_IDs + 1] - offsets[0] and the values array must be at least
   *  this long.
   *
   * \pre n_IDs >= 0
   * \pre values != nullptr
   * \pre offsets != nullptr
   */
  void appendM(const IndexType* values,
               IndexType n_IDs,
               const IndexType* offsets,
               const CellType* AXOM_NOT_USED(types) = nullptr)
  {
    internal::append(n_IDs, values, offsets, m_values, m_offsets);
  }

  /*!
   * \brief Sets the values of the given ID.
   *
   * \param [in] values pointer to the values to set, of length at least
   *  getNumberOfValuesForID( ID ).
   * \param [in] ID the ID of the values to set.
   *
   * \note The number of values associated with the given ID may not be changed,
   *  only the values themselves.
   *
   * \pre ID >= 0 && ID < getNumberOfIDs()
   * \pre values != nullptr
   */
  void set(const IndexType* values, IndexType ID) { setM(values, ID, 1); }

  /*!
   * \brief Sets the values of multiple IDs starting with the given ID.
   *
   * \param [in] values pointer to the values to set, of length at least
   *  the sum of the number of values of each ID to set.
   * \param [in] start_ID the ID to start at.
   * \param [in] n_IDs the number of IDs to set.
   *
   * \note The number of values associated with the given IDs may not be
   *  changed, only the values themselves.
   *
   * \pre start_ID >= 0 && start_ID + n_IDs < getNumberOfIDs()
   * \pre values != nullptr
   */
  void setM(const IndexType* values, IndexType start_ID, IndexType n_IDs)
  {
    internal::set(start_ID, values, n_IDs, m_values, m_offsets);
  }

  /*!
   * \brief Insert the values of a new ID before the given ID.
   *
   * \param [in] values pointer to the values to set, of length at least
   *  n_values.
   * \param [in] start_ID the ID to start at.
   * \param [in] n_values the number of values the new ID has.
   * \param [in] type not used, does not need to be specified.
   *
   * \pre start_ID >= 0 && start_ID <= getNumberOfIDs()
   * \pre n_IDs >= 0
   * \pre values != nullptr
   * \pre offsets != nullptr
   */
  void insert(const IndexType* values,
              IndexType start_ID,
              IndexType n_values,
              CellType AXOM_NOT_USED(type) = UNDEFINED_CELL)
  {
    IndexType offsets[2];
    offsets[0] = 0;
    offsets[1] = n_values;
    insertM(values, start_ID, 1, offsets);
  }

  /*!
   * \brief Insert the values of new IDs before the given ID.
   *
   * \param [in] values pointer to the values to insert.
   * \param [in] start_ID the ID to insert at.
   * \param [in] n_IDs the number of IDs to insert.
   * \param [in] offsets the offsets array of length at least n_IDs + 1.
   * \param [in] types not used, does not need to be specified.
   *
   * \note The number of values to insert is given by
   *  offsets[n_IDs + 1] - offsets[0] and the values array must be at least
   *  this long.
   *
   * \pre start_ID >= 0 && start_ID <= getNumberOfIDs()
   * \pre n_IDs >= 0
   * \pre values != nullptr
   * \pre offsets != nullptr
   */
  void insertM(const IndexType* values,
               IndexType start_ID,
               IndexType n_IDs,
               const IndexType* offsets,
               const CellType* AXOM_NOT_USED(types) = nullptr)
  {
    internal::insert(start_ID, n_IDs, values, offsets, m_values, m_offsets);
  }

  /// @}

private:
  CellType m_cell_type;
  Array<IndexType>* m_values;
  Array<IndexType>* m_offsets;

  DISABLE_COPY_AND_ASSIGNMENT(ConnectivityArray);
  DISABLE_MOVE_AND_ASSIGNMENT(ConnectivityArray);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MINT_ConnectivityArray_indirection_HPP_ */
