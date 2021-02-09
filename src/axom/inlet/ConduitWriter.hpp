// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file ConduitWriter.hpp
 *
 * \brief This file contains the class definition of the ConduitWriter.
 *******************************************************************************
 */

#ifndef INLET_CONDUITWRITER_HPP
#define INLET_CONDUITWRITER_HPP

#include "axom/sidre.hpp"
#include "axom/inlet/Writer.hpp"

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class ConduitWriter
 *
 * \brief A Writer that leverages Conduit's output capabilities to dump an Inlet
 * hierarchy
 *
 * \see Writer
 *******************************************************************************
 */
class ConduitWriter : public Writer
{
public:
  /*!
  *******************************************************************************
  * \brief A constructor for ConduitWriter.
  * 
  * \param [in] fileName The name of the file the documentation should be written to.
  *******************************************************************************
  */
  ConduitWriter(const std::string& filename, const std::string& protocol);

  void documentTable(const Table& table) override;

  void finalize() override;

  virtual ~ConduitWriter() = default;

private:
  conduit::Node m_root;
  const std::string m_filename;
  const std::string m_protocol;
};

}  // namespace inlet
}  // namespace axom

#endif
