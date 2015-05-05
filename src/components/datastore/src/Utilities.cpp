/*!
 ******************************************************************************
 *
 * \file
 *
 * \brief   Implementation file for utility routines.
 *
 ******************************************************************************
 */


// Associated header file
#include "Utilities.hpp"

#include <cstdlib>
#include <sstream>


namespace asctoolkit {

namespace utilities {

/*
 * Routine that prints message, file, line number to std out.
 */
void printMessage(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   std::cout << "File ``" << filename << "'' at line " << line << std::endl;
   std::cout << "MESSAGE: " << std::endl << message << std::endl;
}

/*
 * Aborts the program after printing an error message with file and
 * line number information.
 */
void processAbort(
   const std::string& message,
   const std::string& filename,
   const int line) 
{
   utilities::printMessage( message, filename, line);
   std::cout << "PROGRAM TERMINATION!!!" << std::endl;
   exit(-1);   
}

/*
 * Prints a warning message with file and line number information,
 * but does not abort the program.
 */
void processWarning(
   const std::string& message,
   const std::string& filename,
   const int line)
{
   utilities::printMessage( message, filename, line);
}

/*
 * Convert int to string.
 */
std::string intToString(int val, int min_width)
{
   int tmp_width = (min_width > 0 ? min_width : 1);
   std::ostringstream oss;
   if (val < 0) {
      oss << '-' << std::setw(tmp_width - 1) << std::setfill('0') << -val;
   } else {
      oss << std::setw(tmp_width) << std::setfill('0') << val;
   }
   oss << std::flush;

   return oss.str();  //returns the string form of the stringstream object
}


}  // ending brace for utilities namespace

}  // ending brace for asctoolkit namespace
