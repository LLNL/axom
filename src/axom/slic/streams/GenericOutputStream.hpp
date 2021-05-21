// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file GenericOutputStream.hpp
 *
 */

#ifndef GENERICOUTPUTSTREAM_HPP_
#define GENERICOUTPUTSTREAM_HPP_

#include "axom/slic/core/LogStream.hpp"

#include "axom/core/Macros.hpp"

// C/C++ includes
#include <iostream>  // for ostream

namespace axom
{
namespace slic
{
/*!
 * \class GenericOutputStream
 *
 * \brief Concrete instance of LogStream, which implements functionality for
 *  logging messages to a C++ ostream object, e.g., std::cout, std::cerr, a
 *  std::ofstream, std::ostringstream, etc.
 *
 * \see LogStream Logger
 */
class GenericOutputStream : public LogStream
{
public:
  /*!
   * \brief Constructs a GenericOutpuStream instance with the given stream.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \pre os != NULL
   */
  GenericOutputStream(std::ostream* os);

  /*!
   * \brief Constructs a GenericOutputStream instance with the given stream and
   *  message formatting.
   * \param [in] os pointer to a user-supplied ostream instance.
   * \param [in] format the format string.
   * \see LogStream::setFormatString for the format string.
   */
  GenericOutputStream(std::ostream* os, const std::string& format);

  /*!
   * \brief Destructor.
   */
  virtual ~GenericOutputStream();

  /// \see LogStream::append
  virtual void append(message::Level msgLevel,
                      const std::string& message,
                      const std::string& tagName,
                      const std::string& fileName,
                      int line,
                      bool filter_duplicates);

private:
  std::ostream* m_stream;

  /*!
   * \brief Default constructor.
   * \note Made private to prevent applications from using it.
   */
  GenericOutputStream() : m_stream(static_cast<std::ostream*>(nullptr)) {};

  DISABLE_COPY_AND_ASSIGNMENT(GenericOutputStream);
  DISABLE_MOVE_AND_ASSIGNMENT(GenericOutputStream);
};

} /* namespace slic */

} /* namespace axom */

#endif /* GENERICOUTPUTSTREAM_HPP_ */
