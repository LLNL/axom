// Copyright (c) 2017-2022, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/slic/streams/GenericOutputStream.hpp"

#include "axom/core/Macros.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

namespace axom
{
namespace slic
{
GenericOutputStream::GenericOutputStream(std::ostream* os) : m_stream(os) { }

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(const std::string& stream)
{
  if(stream == "cout")
  {
    m_stream = &std::cout;
  }
  else if(stream == "cerr")
  {
    m_stream = &std::cerr;
  }
  else
  {
    m_stream = new std::ofstream(stream);
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
  // Fix newline and tab characters if needed
  std::string format_fixed = axom::utilities::string::replaceAllInstances(
    axom::utilities::string::replaceAllInstances(format, "\\n", "\n"),
    "\\t",
    "\t");
  this->setFormatString(format_fixed);
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

//------------------------------------------------------------------------------
void GenericOutputStream::flush() { m_stream->flush(); }

} /* namespace slic */

} /* namespace axom */
