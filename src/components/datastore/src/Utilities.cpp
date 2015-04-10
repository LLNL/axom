/*
 * Utilities.cpp
 */

#include "Utilities.hpp"

#include <cstdlib>


namespace DataStoreNS {

namespace utilities {

/*!
 * Routine that prints message, file, line number to std out.
 */
void PrintMessage(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   std::cout << "File ``" << filename << "'' at line " << line << std::endl;
   std::cout << "MESSAGE: " << std::endl << message << std::endl;
}

/*!
 * Aborts the program after printing an error message with file and
 * line number information.
 */
void Abort(
   const std::string& message,
   const std::string& filename,
   const int line) 
{
   utilities::PrintMessage( message, filename, line);
   std::cout << "PROGRAM TERMINATION!!!" << std::endl;
   exit(-1);   
}

/*!
 * Prints a warning message with file and line number information,
 * but does not abort the program.
 */
void Warning(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   utilities::PrintMessage( message, filename, line);
}


}  // ending brace for utilities namespace

}  // ending brace for datastore namespace
