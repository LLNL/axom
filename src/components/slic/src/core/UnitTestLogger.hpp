/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *  \file UnitTestLogger.hpp
 *
 *  \brief Header file containing definition of UnitTestLogger class.
 *
 */

#ifndef UNITTESTLOGGER_HPP_
#define UNITTESTLOGGER_HPP_

// Other CS Toolkit headers
#include "common/ATKMacros.hpp"    // defines DISABLE_{COPY,MOVE}_AND_ASSIGNMENT

// slic component headers
#include "slic.hpp"
#include "GenericOutputStream.hpp"

namespace axom {
namespace slic {

/*!
 * \class UnitTestLogger
 *
 * \brief UnitTestLogger is a simple wrapper around the initialization and
 * finalize operations of the slic::Logger class for CS Toolkit unit tests.
 *
 * To use, create an instance of of this class before tests are run. This
 * initializes the slic logger. When the object is destroyed (e.g., goes out
 * of scope), the slic logger is finalized. For example, when using gtest,
 * a simple main program can be used in each test source file to do this:
 *
 * \verbatim
 *
 *  int main(int argc, char* argv[])
 *  {
 *     int result = 0;
 *
 *     ::testing::InitGoogleTest(&argc, argv);
 *
 *     // create & initialize test logger, finalized when exiting main scope
 *     UnitTestLogger logger;
 *
 *     result = RUN_ALL_TESTS();
 *
 *     return 0;
 *  }
 *
 * \endverbatim
 */
class UnitTestLogger
{
public:

  /*!
   * \brief Constructor initializes slic loging environment.
   */
  UnitTestLogger()
  {
    initialize();
    setLoggingMsgLevel( message::Debug );

    // Formatting for warning, errors and fatal message
    std::string wefFormatStr =
      std:: string( "\n***********************************\n")+
      std:: string( "[<LEVEL> in line <LINE> of file <FILE>]\n") +
      std:: string( "MESSAGE=<MESSAGE>\n" ) +
      std:: string( "***********************************\n");

    // Simple formatting for debug and info messages
    std::string diFormatStr = "[<LEVEL>] <MESSAGE> \n";

    GenericOutputStream* wefStream
      = new GenericOutputStream(&std::cout, wefFormatStr);
    GenericOutputStream* diStream
      = new GenericOutputStream(&std::cout, diFormatStr);

    addStreamToMsgLevel(wefStream,  message::Error);
    addStreamToMsgLevel(wefStream,  message::Warning);
    addStreamToMsgLevel(diStream,   message::Info);
    addStreamToMsgLevel(diStream,   message::Debug);

  }

  /*!
   * \brief Destructor finalizes slic loging environment.
   */
  ~UnitTestLogger()
  {
    finalize();
  }

private:
  DISABLE_COPY_AND_ASSIGNMENT(UnitTestLogger);
  DISABLE_MOVE_AND_ASSIGNMENT(UnitTestLogger);

};

} /* end namespace slic */
} /* end namespace axom */

#endif /* UNITTESTLOGGER_HPP_ */
