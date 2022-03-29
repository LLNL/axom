// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
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
GenericOutputStream::GenericOutputStream(const std::string& stream)
{
  if (stream == "cout")
  {
    m_stream = &std::cout;
  }
  else if (stream == "cerr")
  {
    m_stream = &std::cerr;
  }
  else
  {
    std::ofstream ofs;
    ofs.open(stream);
    m_stream = &ofs;
  }
}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os,
                                         const std::string& format)
  : m_stream(os)
{
  this->setFormatString(format);
}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(const std::string& stream,
                                         const std::string& format)
  : GenericOutputStream::GenericOutputStream(stream)
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
                                 bool AXOM_UNUSED_PARAM(filtered_duplicates))
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
