//Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
//other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
//SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MAP_HPP_
#define AXOM_MAP_HPP_

namespace axom
{


template <typename KType, typename VType>
class Map
{
public:

  AXOM_HOST_DEVICE
  Map() {}

  AXOM_HOST_DEVICE
  VType& operator[](KType query){}

  AXOM_HOST_DEVICE
  const VType& operator[](KType query) const {}

  AXOM_HOST_DEVICE
  Map(Array&& other);

  AXOM_HOST_DEVICE
  Map(Array& other);

};

} /* namespace axom */

#endif /* AXOM_MAP_HPP_ */
