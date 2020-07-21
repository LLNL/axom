// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 *******************************************************************************
 * \file SphinxDocWriter.hpp
 *
 * \brief This file contains the class definition of the SphinxDocWriter.
 *******************************************************************************
 */

#ifndef INLET_SPHINXDOCWRITER_HPP
#define INLET_SPHINXDOCWRITER_HPP

#include <string>
#include <vector>
#include <fstream>

#include "axom/sidre.hpp"
#include "axom/inlet/DocWriter.hpp"

namespace axom
{
namespace inlet
{ 

/*!
 *******************************************************************************
 * \class SphinxDocWriter
 *
 * \brief A DocWriter that is able write documentation in Sphinx RST format for 
 * a given input deck.
 *
 * \see DocWriter
 *******************************************************************************
 */
class SphinxDocWriter : public DocWriter {
public: 
  /*!
  *******************************************************************************
  * \brief A constructor for SphinxDocWriter.
  * 
  * \param [in] fileName The name of the file the documentation should be written to.
  * 
  * \param [in] sidreRootGroup The root of the sidre group that will be traversed
  * to create documentation.
  *
  *******************************************************************************
  */
  SphinxDocWriter(const std::string& fileName, axom::sidre::Group* sidreRootGroup);

  /*!
  *******************************************************************************
  * \brief Writes all documentation to a file.
  * 
  * This generates all RST-syntax documentation and writes it to the pre-specified
  * file.
  *
  *******************************************************************************
  */
  void writeDocumentation();

private:
  /*!
   *****************************************************************************
   * \brief Accumulate the rstTables vector with documentation data.
   *
   * This recursively identifies all Fields and Tables for the documentation.
   *
   * \param [in] sidreGroup The root of the sidre tree to traverse
   *
   *****************************************************************************
   */
  void writeDocumentationHelper(axom::sidre::Group* sidreGroup);

  /*!
   *****************************************************************************
   * \brief Writes the title in RST syntax.
   *
   * This writes a title to the ostringstream in RST syntax.
   *
   * \param [in] title The title to be written
   *
   *****************************************************************************
   */
  void writeTitle(const std::string& title);

  /*!
   *****************************************************************************
   * \brief Writes the sub-title in RST syntax.
   *
   * This writes a sub-title to the ostringstream in RST syntax.
   *
   * \param [in] sub The sub-title to be written
   *
   *****************************************************************************
   */
  void writeSubtitle(const std::string& sub);

  /*!
   *****************************************************************************
   * \brief Writes a 4 column table in RST syntax.
   *
   * This writes a 4 column table to the ostringstream in RST syntax. The number
   * of rows are determined by the number of Fields found in the vector.
   *
   * \param [in] sub The title of the table written
   * 
   * \param [in] rstTable The 2 dimensional vector containing information to
   * be translated into an RST table
   *
   *****************************************************************************
   */
  void writeTable(const std::string& title, 
                  const std::vector<std::vector<std::string>>& rstTable);

  /*!
   *****************************************************************************
   * \brief Writes all tables and their respective titles and descriptions.
   *
   * This parses all of the information from m_rstTables into RST-syntax 
   * documentation and writes it to the ostringstream.
   *
   *****************************************************************************
   */
  void writeAllTables();

  /*!
  *******************************************************************************
  * \struct TableData
  *
  * \brief A struct to store data associated with each inlet::Table.
  *
  *******************************************************************************
  */
  struct TableData {
    /*!
    *******************************************************************************
    * \brief A constructor for the TableData struct
    * 
    * This initializes the RST table's column labels.
    *
    *******************************************************************************
    */
    TableData() {
      rstTable = {{"Field Name", "Description", "Default Value", 
                             "Range/Valid Values", "Required"}};
    }

    std::string tableName;
    std::string description;
    std::vector<std::vector<std::string>> rstTable;
  };

  /*!
    *******************************************************************************
    * \brief Stores all field info into m_rstTables
    * 
    * This extracts information about description, required, default values, 
    * and range from the Sidre Group corresponding to a Field. Then, the information
    * is stored in m_rstTables.
    *
    *******************************************************************************
    */
  void collectFieldInfo(axom::sidre::Group* sidreGroup);

  axom::sidre::Group* m_sidreRootGroup;
  std::ofstream m_outFile;
  std::ostringstream m_oss;
  std::vector<TableData> m_rstTables;
  std::string m_fileName;
  TableData m_currentTable;
};

}
}

#endif
