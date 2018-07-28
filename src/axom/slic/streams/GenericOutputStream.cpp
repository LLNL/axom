/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#include "axom/slic/streams/GenericOutputStream.hpp"

#include "axom/core/Macros.hpp"

namespace axom
{
namespace slic
{

GenericOutputStream::GenericOutputStream( std::ostream* os ) :
  m_stream( os )
{}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os,
                                         const std::string& format) :
  m_stream( os )
{
  this->setFormatString( format );
}

//------------------------------------------------------------------------------
GenericOutputStream::~GenericOutputStream()
{}

//------------------------------------------------------------------------------
void GenericOutputStream::append( message::Level msgLevel,
                                  const std::string& message,
                                  const std::string& tagName,
                                  const std::string& fileName,
                                  int line,
                                  bool AXOM_NOT_USED(filtered_duplicates) )
{
  if ( m_stream == nullptr )
  {
    std::cerr << "ERROR: NULL stream!\n";
    return;
  }

  (*m_stream) << this->getFormatedMessage( message::getLevelAsString(msgLevel),
                                           message,
                                           tagName,
                                           "",
                                           fileName,
                                           line );
}

} /* namespace slic */

} /* namespace axom */
