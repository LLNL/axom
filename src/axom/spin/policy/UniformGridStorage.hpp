// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_POLICY_UGRID_STORAGE_HPP
#define AXOM_SPIN_POLICY_UGRID_STORAGE_HPP

#include "axom/core/Array.hpp"

namespace axom
{
namespace spin
{
namespace policy
{
template <typename T>
struct DynamicGridStorage;

template <typename T>
struct DynamicGridView;

/*!
 * \brief Policy for storing uniform grid as an array of arrays. Supports
 *  efficient single-element insertions and bin clears.
 */
template <typename T>
struct DynamicGridStorage
{
  using ViewType = DynamicGridView<T>;
  using ConstViewType = DynamicGridView<const T>;

  using BinType = axom::Array<T>;
  using ConstBinType = const axom::Array<T>;

  using BinRef = axom::ArrayView<T>;
  using ConstBinRef = axom::ArrayView<const T>;

  DynamicGridStorage(int allocID)
    : m_bins(0, 0, allocID)
    , m_allocatorID(allocID)
  { }

  int getAllocatorID() const { return m_allocatorID; }

  void setNumBins(IndexType nbins) { m_bins.resize(nbins); }
  IndexType getNumBins() const { return m_bins.size(); }

  bool isValidIndex(IndexType index) const
  {
    return index >= 0 && index < getNumBins();
  }

  void initialize(axom::ArrayView<const IndexType> binSizes)
  {
    axom::Array<BinType> host_bins;
    for(int i = 0; i < binSizes.size(); i++)
    {
      host_bins.emplace_back(ArrayOptions::Uninitialized {},
                             binSizes[i],
                             binSizes[i],
                             m_allocatorID);
    }
    m_bins = axom::Array<BinType>(std::move(host_bins), m_allocatorID);
  };

  void insert(IndexType gridIdx, const T& elem)
  {
    m_bins[gridIdx].push_back(elem);
  }

  void clear(IndexType gridIdx)
  {
    if(isValidIndex(gridIdx))
    {
      m_bins[gridIdx].clear();
    }
  }

  // getters
  BinRef getBinContents(IndexType gridIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }
  ConstBinRef getBinContents(IndexType gridIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }

  T& get(IndexType gridIdx, IndexType binIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  const T& get(IndexType gridIdx, IndexType binIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  axom::Array<axom::Array<T>> m_bins;
  int m_allocatorID;
};

template <typename T>
struct DynamicGridView
{
  using BaseT = typename std::remove_const<T>::type;
  using BinType = axom::ArrayView<T>;
  using ConstBinType = axom::ArrayView<T>;

  using BinStoreType = typename std::conditional<std::is_const<T>::value,
                                                 const axom::Array<T>,
                                                 axom::Array<T>>::type;

  DynamicGridView(DynamicGridStorage<BaseT>& in) : m_bins(in.m_bins) { }

  DynamicGridView(const DynamicGridStorage<BaseT>& in) : m_bins(in.m_bins) { }

  AXOM_HOST_DEVICE IndexType getNumBins() const { return m_bins.size(); }

  AXOM_HOST_DEVICE bool isValidIndex(IndexType index) const
  {
    return index >= 0 && index < getNumBins();
  }
  // getters
  AXOM_HOST_DEVICE BinType getBinContents(IndexType gridIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }
  AXOM_HOST_DEVICE ConstBinType getBinContents(IndexType gridIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx];
  }

  T& get(IndexType gridIdx, IndexType binIdx)
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  const T& get(IndexType gridIdx, IndexType binIdx) const
  {
    SLIC_ASSERT(isValidIndex(gridIdx));
    return m_bins[gridIdx][binIdx];
  }

  axom::ArrayView<BinStoreType> m_bins;
};

}  // namespace policy
}  // namespace spin
}  // namespace axom
#endif  // AXOM_SPIN_POLICY_UGRID_STORAGE_HPP
