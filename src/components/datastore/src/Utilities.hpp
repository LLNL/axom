/*
 * Utilities.hpp
 */

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <iostream>
#include <iomanip>

namespace DataStoreNS
{

//-----------------------------------------------------------------------------
//
/// The ASCTK_ERROR macro can be used to capture errors and output useful 
/// messages in the datastore. It currently prints the message, file, and line
/// number and then ends the program by calling exit(). 
/// If we need to call abort() instead, or use this in MPI programs, we can 
/// augment the Error() method  to do the right thing.
///
//-----------------------------------------------------------------------------
#define ASCTK_ERROR( msg )                                      \
do {                                                                \
    std::ostringstream oss;                                         \
    oss << "Error Message: " << msg << std::ends;                   \
    DataStoreNS::utilities::Abort( oss.str(), __FILE__, __LINE__);  \
} while (0)


//-----------------------------------------------------------------------------
//
/// The ASCTK_WARNING macro can be used to output warning messages
/// in the datastore. It the message, file, and line number but does not 
/// force the program to exit.
///
//-----------------------------------------------------------------------------
#define ASCTK_WARNING( msg )                                    \
do {                                                                \
    std::ostringstream oss;                                         \
    oss << "Warning Message: " << msg << std::ends;                 \
    DataStoreNS::utilities::Warning( oss.str(), __FILE__, __LINE__);\
} while (0)


//-----------------------------------------------------------------------------
//
/// The ASCTK_ASSERT macro can be used to capture an assertion when
/// the given expression does not evaluate to true.  After printing
/// the failed assertion, it forces the program to exit in the same way
/// as ASCTK_ERROR above.
///
//-----------------------------------------------------------------------------
#define ASCTK_ASSERT( EXP )                                        \
do {                                                                   \
    if (!(EXP)) {                                                      \
       std::ostringstream oss;                                         \
       oss << "Failed Assert: " << # EXP << std::ends;                 \
       DataStoreNS::utilities::Abort( oss.str(), __FILE__, __LINE__);  \
    }                                                                  \
} while (0)


//-----------------------------------------------------------------------------
//
/// The ASCTK_ASSERT_MSG macro can be used to capture an assertion when
/// the given expression does not evaluate to true.  After printing
/// the failed assertion and associated message, it forces the program to 
/// exit in the same way as ASCTK_ERROR above.
///
//-----------------------------------------------------------------------------
#define ASCTK_ASSERT_MSG( EXP, msg )                               \
do {                                                                   \
    if (!(EXP)) {                                                      \
       std::ostringstream oss;                                         \
       oss << "Failed Assert: " << # EXP << std::endl << msg << std::ends;\
       DataStoreNS::utilities::Abort( oss.str(), __FILE__, __LINE__);  \
    }                                                                  \
} while (0)


namespace utilities {

   /*!
    * Routine that prints message, file, line number to std out.
    */
   void PrintMessage(
      const std::string& message,
      const std::string& filename,
      const int line);

   /*!
    * Aborts the program after printing an error message with file and
    * line number information.
    */
   void Abort(
      const std::string& message,
      const std::string& filename,
      const int line); 

   /*!
    * Prints a warning message with file and line number information,
    * but does not abort the program.
    */
   void Warning(
      const std::string& message,
      const std::string& filename,
      const int line);

}  // ending brace for utilities namespace


}  // ending brace for datastore namespace

#endif /* UTILITIES_HPP_ */
