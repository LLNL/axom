// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ConduitWriter.cpp
 *
 * \brief This file contains the class implementation of the ConduitWriter.
 *******************************************************************************
 */

#include "axom/inlet/ConduitWriter.hpp"

#include <fstream>

#include "axom/core/utilities/FileUtilities.hpp"
#include "axom/core/utilities/StringUtilities.hpp"
#include "axom/inlet/inlet_utils.hpp"
#include "axom/inlet/Table.hpp"

#include "fmt/fmt.hpp"
#include "axom/slic.hpp"

namespace axom
{
namespace inlet
{
ConduitWriter::ConduitWriter(const std::string& filename, const std::string& protocol) : m_filename(filename), m_protocol(protocol)
{
}

void ConduitWriter::documentTable(const Table& table)
{
  // conduit::Node* node = &m_root;
  // std::vector<std::string> tokens;
  // axom::utilities::string::split(tokens, table.name(), '/');
  // for(const auto& token : tokens)
  // {
  //   int token_as_int;
  //   bool is_int = checkedConvertToInt(token, token_as_int);
  //   if(is_int && token_as_int < node->number_of_children())
  //   {
  //     node = &((*node)[token_as_int]);
  //   }
  //   else
  //   {
  //     node = &((*node)[token]);
  //   }
  // }

  // conduit::Node& table_node = *node;
  for(const auto& field_entry : table.getChildFields())
  {
    if (*(field_entry.second))
    {
      conduit::Node& field_node = m_root[field_entry.first];
      const auto sidreGroup = field_entry.second->sidreGroup();
      sidreGroup->createNativeLayout(field_node);
    }
  }
}

void ConduitWriter::finalize()
{
  SLIC_INFO(m_root.to_yaml());
}



}  // end namespace inlet
}  // end namespace axom
