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
 * \file LogStream.cpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "LogStream.hpp"

#include <sstream>

namespace asctoolkit {

namespace logapi {


//------------------------------------------------------------------------------
LogStream::LogStream() :
    m_formatString("*****\n[<LEVEL>]\n\n <MESSAGE> \n\n <FILE>\n<LINE>\n****\n")
{

}

//------------------------------------------------------------------------------
LogStream::~LogStream()
{

}

//------------------------------------------------------------------------------
void LogStream::replaceKey( std::string& msg,
                            const std::string& key,
                            const std::string& value )
{

  std::size_t pos = msg.find( key );

  if ( pos != std::string::npos ) {

    msg = msg.substr(0,pos) +
          value +
          msg.substr(pos+key.length(),msg.length()-1);

  } // END if

}

//------------------------------------------------------------------------------
std::string LogStream::getFormatedMessage( const std::string& msgLevel,
                                           const std::string& message,
                                           const std::string& tagName,
                                           const std::string& fileName,
                                           int line )
{
  std::string msg = m_formatString;

  this->replaceKey( msg, "<LEVEL>", msgLevel );
  this->replaceKey( msg, "<MESSAGE>", message );
  this->replaceKey( msg, "<TAG>", tagName );
  this->replaceKey( msg, "<FILE>", fileName );

  if ( line != MSG_IGNORE_LINE ) {

    std::ostringstream oss;
    oss << line;

    this->replaceKey( msg, "<LINE>", oss.str() );

  } else {

    this->replaceKey( msg, "<LINE>", "" );

  }

  return( msg );
}

} /* namespace logapi */

} /* namespace asctoolkit */
