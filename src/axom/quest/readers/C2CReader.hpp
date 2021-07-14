// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_C2CREADER_HPP_
#define QUEST_C2CREADER_HPP_

#include "axom/config.hpp"

// This file is only applicable if axom was configured with the C2C library
#ifdef AXOM_USE_C2C

  #include "axom/mint.hpp"
  #include "c2c/C2C.hpp"

  #include <string>
  #include <vector>

namespace axom
{
namespace quest
{
/*
 * \class C2CReader
 *
 * \brief A class to help with reading C2C contour files and converting them to 1D meshes
 *
 * We treat all contours as NURBS curves and sample them to create segment meshes.
 */
class C2CReader
{
public:
  C2CReader() = default;

  virtual ~C2CReader() = default;

  void setFileName(const std::string& fileName) { m_fileName = fileName; }

  void setLengthUnit(c2c::LengthUnit lengthUnit) { m_lengthUnit = lengthUnit; }

  /// Clears data associated with this reader
  void clear();

  /*!
   * \brief Read the contour file provided by \a setFileName()
   * 
   * \return 0 for a successful read; non-zero otherwise
   */
  virtual int read();

  /// \brief Utility function to log details about the read in file
  virtual void log();

  /*!
   * \brief Projects high-order NURBS contours onto a linear mesh using \a segmentsPerPiece linear segments
   * per \a Piece of the contour
   */
  void getLinearMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh,
                     int segmentsPerPiece);

protected:
  int readContour();

protected:
  std::string m_fileName;
  c2c::LengthUnit m_lengthUnit {c2c::LengthUnit::cm};

  std::vector<c2c::NURBSData> m_nurbsData;
};

}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_C2C

#endif  // QUEST_C2CREADER_HPP_
