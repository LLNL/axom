// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file MapBackend.hpp
 *
 * \brief This file contains the class definition of MapBackend.
 *******************************************************************************
 */

#ifndef INLET_BACKEND_HPP
#define INLET_BACKEND_HPP

#include <string>

#include "axom/inlet/Field.hpp"

namespace axom
{
namespace inlet
{

class Backend
{
public:
  virtual void add(Field* field) = 0;
  virtual Field* get(const std::string& name) = 0;
};

} // end namespace inlet
} // end namespace axom

#endif
