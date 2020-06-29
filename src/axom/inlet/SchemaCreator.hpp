// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef INLET_SCHEMACREATOR_HPP
#define INLET_SCHEMACREATOR_HPP

#include <memory>
#include <string>

#include "axom/inlet/Field.hpp"

namespace axom
{
namespace inlet
{

class Table;

class SchemaCreator
{
public:
  // Functions that define the input deck schema
  virtual std::shared_ptr<Table> addTable(const std::string& name,
                                          const std::string& description) = 0;

  virtual std::shared_ptr<Field> addBool(const std::string& name,
                                         const std::string& description) = 0;
  virtual std::shared_ptr<Field> addDouble(const std::string& name,
                                           const std::string& description) = 0;
  virtual std::shared_ptr<Field> addInt(const std::string& name,
                                        const std::string& description) = 0;
  virtual std::shared_ptr<Field> addString(const std::string& name,
                                           const std::string& description) = 0;
};

} // end namespace inlet
} // end namespace axom

#endif
