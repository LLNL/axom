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

#ifndef SLIM_BACKEND_HPP
#define SLIM_BACKEND_HPP

#include <string>

#include "axom/slim/Field.hpp"

namespace axom
{
namespace slim
{

class Backend
{
public:
    virtual void add(Field* field) = 0;
    virtual Field* get(std::string name) = 0;
};

} // end namespace slim
} // end namespace axom

#endif
