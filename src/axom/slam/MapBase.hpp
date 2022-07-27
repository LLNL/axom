// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/**
 * \file MapBase.hpp
 *
 * \brief Contains an Abstract class MapBase
 *
 */

#ifndef SLAM_MAPBASE_HPP_
#define SLAM_MAPBASE_HPP_

#include "axom/core/Macros.hpp"
#include "axom/core/Types.hpp"

#include "axom/slam/Set.hpp"

namespace axom
{
namespace slam
{
/**
 * \class   MapBase
 *
 * \brief   A base class specifying methods that support operations of a Map,
 *          associating value(s) to each element in a Set. MapBase can be used
 *          as a base class pointer to a templated Map object.
 * \see     Map
 *
 */

template <typename SetPositionType = slam::DefaultPositionType>
class MapBase
{
public:
  using SetPosition = SetPositionType;

public:
  virtual ~MapBase() {};

  /**
   * \brief Get the number of entities in the set used by this map
   * \return The number of entities in the set used in the map.
   */
  virtual SetPosition size() const = 0;

  /**
   * \brief Checks whether the Map is valid.
   * \return   True if valid, false otherwise.
   */
  virtual bool isValid(bool verboseOutput) const = 0;
};

template <typename MapType>
class MapVirtualProxy : public MapBase<typename MapType::SetPosition>
{
public:
  using SetPosition = typename MapType::SetPosition;

public:
  MapVirtualProxy(MapType inst) : m_impl(std::move(inst)) { }

  SetPosition size() const override { return m_impl.size(); }

  bool isValid(bool verboseOutput) const override
  {
    return m_impl.isValid(verboseOutput);
  }

  const MapType& get() const { return m_impl; }
  MapType& get() { return m_impl; }

private:
  MapType m_impl;
};

template <typename DerivedMap, typename SetPositionType = typename DerivedMap::SetPosition>
std::unique_ptr<MapBase<SetPositionType>> makeVirtualMap(DerivedMap value)
{
  auto* vptr = new MapVirtualProxy<DerivedMap>(std::move(value));
  std::unique_ptr<MapBase<SetPositionType>> ret(vptr);
  return ret;
}

}  // end namespace slam
}  // end namespace axom

#endif  // SLAM_MAPBASE_HPP_
