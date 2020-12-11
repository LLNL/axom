// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file Reader.hpp
 *
 * \brief This file contains the pure virtual base class definition of Reader.
 *******************************************************************************
 */

#ifndef INLET_READER_HPP
#define INLET_READER_HPP

#include <string>
#include <unordered_map>
#include <vector>

namespace axom
{
namespace inlet
{
/*!
 *******************************************************************************
 * \class Reader
 *
 * \brief Abstract base class defining the interface of all Reader
 *  classes.
 *
 *  Concrete instances need to inherit from this class and implement these
 *  functions.You will need to add your Reader using Reader::initialize
 *
 * \see LuaReader
 *******************************************************************************
 */
class Reader
{
public:
  /*!
   *****************************************************************************
   * \brief Virtual destructor.
   *****************************************************************************
   */
  virtual ~Reader() {};

  /*!
   *****************************************************************************
   * \brief Parses the given input file.
   *
   * This performs any setup work and parses the given input file.
   * It is required that this is called before using the Reader and overrides
   * any state that was previously there.
   *
   * \param [in] filePath The Input file to be read
   *
   * \return true if the input file was able to be parsed
   *****************************************************************************
   */
  virtual bool parseFile(const std::string& filePath) = 0;

  /*!
   *****************************************************************************
   * \brief Parses the given string.
   *
   * This performs any setup work and parses the given string.
   * It is required that this is called before using the Reader and overrides
   * any state that was previously there.
   *
   * \param [in] inputString The Input file to be read
   *
   * \return true if the string was able to be parsed
   *****************************************************************************
   */
  virtual bool parseString(const std::string& inputString) = 0;

  /*!
   *****************************************************************************
   * \brief Return a boolean out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in] id The identifier to the bool that will be retrieved
   * \param [out] value The value of the bool that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getBool(const std::string& id, bool& value) = 0;

  /*!
   *****************************************************************************
   * \brief Return a double out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the double that will be retrieved
   * \param [out] value The value of the double that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getDouble(const std::string& id, double& value) = 0;

  /*!
   *****************************************************************************
   * \brief Return a int out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the int that will be retrieved
   * \param [out] value The value of the int that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getInt(const std::string& id, int& value) = 0;

  /*!
   *****************************************************************************
   * \brief Return a string out of the input file
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] value The value of the string that was retrieved
   *
   * \return true if the variable was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getString(const std::string& id, std::string& value) = 0;

  /*!
   *****************************************************************************
   * \brief Get an index-integer mapping for the given array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the ints that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getIntMap(const std::string& id,
                         std::unordered_map<int, int>& values) = 0;
  /// \overload
  virtual bool getIntMap(const std::string& id,
                         std::unordered_map<std::string, int>& values) = 0;

  /*!
   *****************************************************************************
   * \brief Get an index-bool mapping for the given array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the bools that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getBoolMap(const std::string& id,
                          std::unordered_map<int, bool>& values) = 0;
  /// \overload
  virtual bool getBoolMap(const std::string& id,
                          std::unordered_map<std::string, bool>& values) = 0;

  /*!
   *****************************************************************************
   * \brief Get an index-double mapping for the given array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the doubles that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getDoubleMap(const std::string& id,
                            std::unordered_map<int, double>& values) = 0;
  /// \overload
  virtual bool getDoubleMap(const std::string& id,
                            std::unordered_map<std::string, double>& values) = 0;

  /*!
   *****************************************************************************
   * \brief Get an index-string mapping for the given Lua array
   *
   * This performs any necessary retrieval and mapping from the given identifier
   * to what is in the input file.
   *
   * \param [in]  id    The identifier to the string that will be retrieved
   * \param [out] map The values of the strings that were retrieved
   *
   * \return true if the array was able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getStringMap(const std::string& id,
                            std::unordered_map<int, std::string>& values) = 0;
  /// \overload
  virtual bool getStringMap(const std::string& id,
                            std::unordered_map<std::string, std::string>& values) = 0;

  /*!
   *****************************************************************************
   * \brief Get the list of indices for an container
   *
   * \param [in]  id    The identifier to the container that will be retrieved
   * \param [out] indices The values of the indices that were retrieved
   *
   * \return true if the indices were able to be retrieved from the file
   *****************************************************************************
   */
  virtual bool getIndices(const std::string& id, std::vector<int>& indices) = 0;
  /// \overload
  virtual bool getIndices(const std::string& id,
                          std::vector<std::string>& indices) = 0;
};

}  // end namespace inlet
}  // end namespace axom

#endif
