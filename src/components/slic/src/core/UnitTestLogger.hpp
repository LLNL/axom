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
 * To use, creste an instance of of this class before tests are run. This
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

     std::string format = 
        std::string("\n***********************************\n")+
        std::string( "LEVEL=<LEVEL>\n" ) +
        std::string( "MESSAGE=<MESSAGE>\n" ) +
        std::string( "FILE=<FILE>\n" ) +
        std::string( "LINE=<LINE>\n" ) +
        std::string("***********************************\n");

     setLoggingMsgLevel( message::Debug );
     addStreamToAllMsgLevels(new GenericOutputStream(&std::cout, format));
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
#ifdef USE_CXX11
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
