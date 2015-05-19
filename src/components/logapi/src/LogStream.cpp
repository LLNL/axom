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
    m_fmtString( "[<MTYPE>] <FILE>:<LINE>\n MESSAGE: <MESSAGE>\n" )
{

}

//------------------------------------------------------------------------------
LogStream::~LogStream()
{

}

//------------------------------------------------------------------------------
std::string LogStream::getFormatedMessage( const std::string& msgTypeName,
                                           const std::string& message,
                                           const std::string& fileName,
                                           int line )
{
  std::string msg = m_fmtString;

  std::size_t pos = msg.find( "<MTYPE>" );
  if ( pos != std::string::npos ) {

    msg = msg.substr(0,pos) + msgTypeName + msg.substr(pos+7,msg.length()-1);

  } // END if

  pos = msg.find( "<FILE>" );
  if ( pos != std::string::npos ) {

    msg = msg.substr(0,pos) + fileName + msg.substr(pos+6,msg.length()-1);

  } // END if

  pos = msg.find( "<LINE>" );
  if ( pos != std::string::npos ) {

    std::ostringstream oss;
    oss << line;

    msg = msg.substr(0,pos) + oss.str() + msg.substr(pos+6,msg.length()-1);

  } // END if

  pos = msg.find( "<MESSAGE>" );
  if ( pos != std::string::npos ) {

    msg = msg.substr(0,pos) + message + msg.substr(pos+9,msg.length()-1);

  } // END if

  return( msg );
}

} /* namespace logapi */

} /* namespace asctoolkit */
