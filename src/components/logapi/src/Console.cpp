/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 */

/*!
 *******************************************************************************
 * \file Console.cpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#include "Console.hpp"

// C/C++ includes
#include <iostream>

namespace asctoolkit {

namespace logapi {

Console::Console()
{

}

//------------------------------------------------------------------------------
Console::~Console()
{

}

//------------------------------------------------------------------------------
void Console::append( MessageType msgType,
                      const std::string& msgTypeName,
                      const std::string& message,
                      const std::string& fileName,
                      int line )
{
  std::cout << this->getFormatedMessage( msgTypeName, message, fileName, line );
  std::cout << "\n";
}


} /* namespace logapi */

} /* namespace asctoolkit */
