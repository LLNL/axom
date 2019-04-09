// Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "Relation.hpp"


namespace axom
{
namespace slam
{

/**
 * \brief Definition of static instance of nullSet for all relations
 * \note Should this be a singleton or a global object?  Should the scope be
 *  public?
 */
NullSet Relation::s_nullSet;

} // namespace slam
} // namespace axom
