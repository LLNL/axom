// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/GenericOutputStream.hpp"

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slic
{
GenericOutputStream::GenericOutputStream(std::ostream* os) : m_stream(os) { }

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os,
                                         const std::string& format)
  : m_stream(os)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
GenericOutputStream::~GenericOutputStream() { }

//------------------------------------------------------------------------------
void GenericOutputStream::append(message::Level msgLevel,
                                 const std::string& message,
                                 const std::string& tagName,
                                 const std::string& fileName,
                                 int line,
                                 bool AXOM_NOT_USED(filtered_duplicates))
{
  if(m_stream == nullptr)
  {
    std::cerr << "ERROR: NULL stream!\n";
    return;
  }

  (*m_stream) << this->getFormatedMessage(message::getLevelAsString(msgLevel),
                                          message,
                                          tagName,
                                          "",
                                          fileName,
                                          line);
}

} /* namespace slic */

} /* namespace axom */
