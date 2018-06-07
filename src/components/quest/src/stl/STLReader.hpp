/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef STLREADER_HPP_
#define STLREADER_HPP_

// Axom includes
#include "axom/Macros.hpp"  // for axom macros

// Mint includes
#include "mint/UnstructuredMesh.hpp"

// C/C++ includes
#include <string>  // for std::string
#include <vector>  // for std::vector

namespace axom
{
namespace quest
{

/*!
 * \class STLReader
 *
 * \brief A simple reader for an STL file encoded in the ascii or binary format.
 *
 * STL (STereoLithography) is a common file format for triangle meshes.
 * It encodes a "soup of triangles" by explicitly listing the coordinate
 * positions of the three vertices of each triangle.
 */
class STLReader
{
public:

  /*!
   * \brief Constructor.
   */
  STLReader();

  /*!
   * \brief Destructor.
   */
  virtual ~STLReader();

  /*!
   * \brief Sets the name of the file to read.
   * \param [in] fileName the name of the file to read.
   */
  void setFileName( const std::string& fileName ) { m_fileName = fileName; };

  /*!
   * \brief Clears all internal data-structures
   */
  void clear();

  /*!
   * \brief Reads in the surface mesh from an STL file.
   * \pre m_fileName != ""
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  virtual int read();

  /*!
   * \brief Stores the STL data in the supplied unstructured mesh object.
   * \param [in,out] mesh pointer to the unstructured mesh.
   * \pre mesh != AXOM_NULLPTR.
   */
  void getMesh( mint::UnstructuredMesh< mint::SINGLE_SHAPE >* mesh );

private:
  /*!
   * \brief A predicate to check if the file is in ascii format
   *
   * We can test the size of the STL file to determine if it is in
   * the binary or ascii format.  The binary format is defined
   * to have an 80 byte header followed by a 4 byte integer
   * encoding the number of triangles, followed by the triangle data
   * (50 bytes per triangle).
   *
   * \return True, if the file is ascii encoded, False if it is binary
   */
  bool isAsciiFormat() const;

  /*!
   * \brief Reads an ascii-encoded STL file into memory
   * \note The filename should be set with STLReader::setFileName()
   */
  int readAsciiSTL();

  /*!
   * \brief Reads a binary-encoded STL file into memory
   * \note The filename should be set with STLReader::setFileName()
   */
  int readBinarySTL();

protected:
  std::string m_fileName;

  axom::mint::IndexType m_num_nodes;
  axom::mint::IndexType m_num_faces;

  std::vector<double> m_nodes;

private:

  DISABLE_COPY_AND_ASSIGNMENT(STLReader);
  DISABLE_MOVE_AND_ASSIGNMENT(STLReader);
};

} // end namespace quest
} // end namespace axom

#endif /* STLREADER_HPP_ */
