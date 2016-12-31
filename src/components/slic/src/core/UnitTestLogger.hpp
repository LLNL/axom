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
#include "common/config.hpp"    // defines ATK_USE_CXX11

// slic component headers
#include "slic.hpp"
#include "GenericOutputStream.hpp"

namespace asctoolkit
{
namespace slic
{

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
        std::string("\n***********************************\n")+
        std::string( "[<LEVEL> in line <LINE> of file <FILE>]\n") +
        std::string( "MESSAGE=<MESSAGE>\n" ) +
        std::string("***********************************\n");

     // Simple formatting for debug and info messages
     std::string diFormatStr = "[<LEVEL>] <MESSAGE> \n";

     GenericOutputStream* wefStream
             = new GenericOutputStream(&std::cout, wefFormatStr);
     GenericOutputStream* diStream
             = new GenericOutputStream(&std::cout, diFormatStr);


     addStreamToMsgLevel(wefStream, message::Fatal) ;
     addStreamToMsgLevel(wefStream, message::Error);
     addStreamToMsgLevel(wefStream, message::Warning);
     addStreamToMsgLevel(diStream,  message::Info);
     addStreamToMsgLevel(diStream,  message::Debug);

 }

  /*!
   * \brief Destructor finalizes slic loging environment.
   */
  ~UnitTestLogger()
  {
     finalize();
  }

private:
  //
  // Unimplemented copy ctors and copy-assignment operators.
  //
#ifdef ATK_USE_CXX11
  UnitTestLogger( const UnitTestLogger& source ) = delete;
  UnitTestLogger( UnitTestLogger&& source ) = delete;

  UnitTestLogger& operator=( const UnitTestLogger& rhs ) = delete;
  UnitTestLogger& operator=( const UnitTestLogger&& rhs ) = delete;
#else
  UnitTestLogger( const UnitTestLogger& source );
  UnitTestLogger& operator=( const UnitTestLogger& rhs );
#endif

};


} /* end namespace slic */
} /* end namespace asctoolkit */


#endif /* UNITTESTLOGGER_HPP_ */
