/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file GenericOutputStream.cpp
 *
 * \date May 21, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "GenericOutputStream.hpp"

#include "common/ATKMacros.hpp"

namespace axom {
namespace slic {

GenericOutputStream::GenericOutputStream( std::ostream* os ):
    m_stream( os )
{

}

//------------------------------------------------------------------------------
GenericOutputStream::GenericOutputStream(std::ostream* os, const std::string& format):
    m_stream( os )
{
  this->setFormatString( format );
}

//------------------------------------------------------------------------------
GenericOutputStream::~GenericOutputStream()
{

}

//------------------------------------------------------------------------------
void GenericOutputStream::append( message::Level msgLevel,
                                  const std::string& message,
                                  const std::string& tagName,
                                  const std::string& fileName,
                                  int line,
                                  bool ATK_NOT_USED(filtered_duplicates) )
{
  if ( m_stream == ATK_NULLPTR ) {
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
