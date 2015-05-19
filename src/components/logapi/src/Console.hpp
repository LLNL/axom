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
 * \file Console.hpp
 *
 * \date May 7, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef CONSOLE_HPP_
#define CONSOLE_HPP_

#include "logapi/LogStream.hpp"

namespace asctoolkit {

namespace logapi {

/*!
 *******************************************************************************
 * \class Console
 *
 * \brief A concrete instance of LogStream that implements functionality for
 *  logging to standard out.
 *
 * \see LogStream Logger
 *******************************************************************************
 */
class Console : public LogStream
{
public:
  Console();
  virtual ~Console();

  /// \see LogStream::append
  virtual void append( message::Level msgType,
                       const std::string& message,
                       const std::string& fileName,
                       int line );
private:

  /// \name Disabled Methods
  ///@{

  Console( const Console& ); // Not implemented
  Console& operator=( const Console& ); // Not implemented

  ///@}
};

} /* namespace logapi */
} /* namespace asctoolkit */

#endif /* CONSOLE_HPP_ */
