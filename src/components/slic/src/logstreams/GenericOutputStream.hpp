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
 * \file GenericOutputStream.hpp
 *
 * \date May 21, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef GENERICOUTPUTSTREAM_HPP_
#define GENERICOUTPUTSTREAM_HPP_

#include "slic/LogStream.hpp"

#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp" // for ATK_NULLPTR

// C/C++ includes
#include <iostream> // for ostream

namespace asctoolkit {

namespace slic {

/*!
 *******************************************************************************
 * \class GenericOutputStream
 *
 * \brief Concrete instance of LogStream, which implements functionality for
 *  logging messages to a C++ ostream object, e.g., std::cout, std::cerr, a
 *  std::ofstream, std::ostringstream, etc.
 *
 * \see LogStream Logger
 *******************************************************************************
 */
class GenericOutputStream : public LogStream
{
public:

  /*!
   *****************************************************************************
   * \brief Constructs a GenericOutpuStream instance with the given stream.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \pre os != NULL
   *****************************************************************************
   */
  GenericOutputStream( std::ostream* os );

  /*!
   *****************************************************************************
   * \brief Constructs a GenericOutputStream instance with the given stream and
   *  message formatting.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \param [in] format the format string.
   * \see LogStream::setFormatString for the format string.
   *****************************************************************************
   */
  GenericOutputStream( std::ostream* os, std::string& format );

  /*!
   *****************************************************************************
   * \brief Destructor.
   *****************************************************************************
   */
  virtual ~GenericOutputStream();

  /// \see LogStream::append
  virtual void append( message::Level msgLevel,
                       const std::string& message,
                       const std::string& tagName,
                       const std::string& fileName,
                       int line,
                       bool filter_duplicates );
private:

  std::ostream* m_stream;

  /*!
   *****************************************************************************
   * \brief Default constructor.
   * \note Made private to prevent applications from using it.
   *****************************************************************************
   */
  GenericOutputStream(): m_stream( static_cast< std::ostream* >( ATK_NULLPTR ) ){};

  DISABLE_COPY_AND_ASSIGNMENT(GenericOutputStream);

};

} /* namespace slic */

} /* namespace asctoolkit */

#endif /* GENERICOUTPUTSTREAM_HPP_ */
