// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef CORE_STRING_UTILITIES_H_
#define CORE_STRING_UTILITIES_H_

#include <sstream>
#include <string>
#include <vector>

namespace axom
{
namespace utilities
{
namespace string
{
/*!
 * \brief Tests whether a string ends with another string
 * https://stackoverflow.com/questions/20446201/how-to-check-if-string-ends-with-txt/20446257
 * \param [in] str string to be searched
 * \param [in] suffix string to be checked for
 * \return     boolean value whether suffix was found at the end of str
 */
inline bool endsWith(const std::string& str, const std::string& suffix)
{
  return str.size() >= suffix.size() &&
    0 == str.compare(str.size() - suffix.size(), suffix.size(), suffix);
}

/*!
 * \brief Tests whether a string ends with a char
 * \param [in] str string to be searched
 * \param [in] suffix char to be checked for
 * \return     boolean value whether suffix was found at the end of str
 */
inline bool endsWith(const std::string& str, const char suffix)
{
  return str.size() >= 1 && str[str.size() - 1] == suffix;
}

/*!
 * \brief Tests whether a string starts with another string
 * \param [in] str string to be searched
 * \param [in] prefix string to be checked for
 * \return     boolean value whether prefix was found at start of str
 */
inline bool startsWith(const std::string& str, const std::string& prefix)
{
  return str.size() >= prefix.size() &&
    0 == str.compare(0, prefix.size(), prefix);
}

/*!
 * \brief Tests whether a string starts with a char
 * \param [in] str string to be searched
 * \param [in] prefix char to be checked for
 * \return     boolean value whether prefix was found at start of str
 */
inline bool startsWith(const std::string& str, const char prefix)
{
  return str.size() >= 1 && str[0] == prefix;
}

/*!
 * \brief Splits the given string based on the given delimiter
 * \param [out] tokens    vector that the found tokens are appended to
 * \param [in]  str       string to be tokenized
 * \param [in]  delimiter char to split string on
 */
void split(std::vector<std::string>& tokens,
           const std::string& str,
           const char delimiter);

/*!
 * \brief Converts a string to lowercase
 * \param [inout] str    string to be converted
 */
void toLower(std::string& str);

/*!
 * \brief Converts a string to uppercase
 * \param [inout] str    string to be converted
 */
void toUpper(std::string& str);

/**
 * \brief Retrieves the last "n" tokens of a string split with the specified delimiter
 * \param [in] input The string to split
 * \param [in] n The number of tokens to retrieve
 * \param [in] delim The delimiter to split with
 * 
 * \return A list of tokens (of size \a n )
 * 
 * Splits a string starting from the end of the string into a maximum of \a n tokens
 */
std::vector<std::string> splitLastNTokens(const std::string& input,
                                          const std::size_t n,
                                          const char delim);

}  // end namespace string
}  // end namespace utilities
}  // end namespace axom

#endif  //  CORE_STRING_UTILITIES_H_
