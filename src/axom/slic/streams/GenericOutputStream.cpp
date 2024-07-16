// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
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
GenericOutputStream::GenericOutputStream(std::ostream* os)
  : m_stream(os)
  , m_file_name()
  , m_opened(true)
{ }

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(const std::string& stream)
{
  if(stream == "cout")
  {
    m_stream = &std::cout;
    m_file_name = std::string();
    m_opened = true;
  }
  else if(stream == "cerr")
  {
    m_stream = &std::cerr;
    m_file_name = std::string();
    m_opened = true;
  }
  else
  {
    m_stream = new std::ofstream();
    m_file_name = stream;
    m_opened = false;
  }
}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os,
                                         const std::string& format)
  : m_stream(os)
  , m_file_name()
  , m_opened(true)
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
                                 bool AXOM_UNUSED_PARAM(filtered_duplicates),
                                 bool AXOM_UNUSED_PARAM(tag_stream_only))
{
  if(m_stream == nullptr)
  {
    std::cerr << "ERROR: NULL stream!\n";
    return;
  }

  if(!m_opened)
  {
    std::ofstream* ofs = dynamic_cast<std::ofstream*>(m_stream);
    if(ofs != nullptr)
    {
      ofs->open(m_file_name);
      m_opened = true;
    }
  }

  (*m_stream) << this->getFormatedMessage(message::getLevelAsString(msgLevel),
                                          message,
                                          tagName,
                                          "",
                                          "",
                                          fileName,
                                          line);
}

//------------------------------------------------------------------------------
void GenericOutputStream::outputLocal() { m_stream->flush(); }

//------------------------------------------------------------------------------
void GenericOutputStream::flush() { m_stream->flush(); }

} /* namespace slic */

} /* namespace axom */
