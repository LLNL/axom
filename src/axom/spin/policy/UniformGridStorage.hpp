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

/*!
 * \brief Policy for storing uniform grid as an array of arrays. Supports
 *  efficient single-element insertions and bin clears.
 */
template <typename T>
struct DynamicGridStorage
{
  using BinType = axom::ArrayView<T>;
  using ConstBinType = axom::ArrayView<const T>;
  void setNumBins(IndexType nbins) { m_bins.resize(nbins); }
  IndexType getNumBins() const { return m_bins.size(); }

  void initialize(axom::ArrayView<const IndexType> binSizes)
  {
    for(int i = 0; i < binSizes.size(); i++)
    {
      m_bins[i].resize(binSizes[i]);
    }
  };

  void insert(IndexType gridIdx, const T& elem)
  {
    m_bins[gridIdx].push_back(elem);
  }

  void clear(IndexType gridIdx) { m_bins[gridIdx].clear(); }

  // getters
  BinType getBin(IndexType gridIdx) { return m_bins[gridIdx]; }
  ConstBinType getBin(IndexType gridIdx) const { return m_bins[gridIdx]; }

  T& get(IndexType gridIdx, IndexType binIdx)
  {
    return m_bins[gridIdx][binIdx];
  }

  const T& get(IndexType gridIdx, IndexType binIdx) const
  {
    return m_bins[gridIdx][binIdx];
  }

  axom::Array<axom::Array<T>> m_bins;
};

}  // namespace policy
}  // namespace spin
}  // namespace axom
#endif  // AXOM_SPIN_POLICY_UGRID_STORAGE_HPP
