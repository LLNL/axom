// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_FIELD_HPP
#define INLET_FIELD_HPP

#include "axom/sidre.hpp"

#include <memory>

namespace axom
{
namespace inlet
{

class Field : public std::enable_shared_from_this<Field>
{
public:
  Field(axom::sidre::Group* sidreGroup) :
    m_sidreGroup(sidreGroup) {}
  axom::sidre::Group* sidreGroup() { return m_sidreGroup; };

  std::shared_ptr<Field> required(bool isRequired);
  bool required();
private:
  // This Field's sidre group
  axom::sidre::Group* m_sidreGroup = nullptr;
};

} // end namespace inlet
} // end namespace axom

#endif
