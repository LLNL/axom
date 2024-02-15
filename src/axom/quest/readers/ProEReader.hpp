// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef QUEST_PROEREADER_HPP_
#define QUEST_PROEREADER_HPP_

// Axom includes
#include "axom/config.hpp"
#include "axom/core/Macros.hpp"
#include "axom/mint/mesh/UnstructuredMesh.hpp"
#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/slam/BitSet.hpp"

// C/C++ includes
#include <string>  // for std::string
#include <vector>  // for std::vector

namespace axom
{
namespace quest
{
/*!
 * \class ProEReader
 *
 * \brief A simple reader for a Pro/E file encoded in ascii.
 *        Current support for only tetrahedra.
 *
 * Pro/Engineer (also known as Creo) is a modeling application.
 *
 * \note Pro/E node IDs start at 1. ProEReader adjusts and stores
 *  node IDs to start at 0 for indexing.
 */
class ProEReader
{
public:
  constexpr static int NUM_NODES_PER_TET = 4;
  constexpr static int NUM_COMPS_PER_NODE = 3;
  using Point3D = primal::Point<double, NUM_COMPS_PER_NODE>;
  using BBox3D = primal::BoundingBox<double, 3>;
  using TetPred = std::function<bool(int[4], int, std::vector<double>&)>;

public:
  /*!
   * \brief Constructor.
   */
  ProEReader();

  /*!
   * \brief Destructor.
   */
  virtual ~ProEReader();

  /*!
   * \brief Sets the name of the file to read.
   * \param [in] fileName the name of the file to read.
   */
  void setFileName(const std::string& fileName) { m_fileName = fileName; };

  /*!
   * \brief Returns the number of nodes of the mesh.
   * \return numNodes the number of nodes.
   */
  int getNumNodes() const { return m_num_nodes; };

  /*!
   * \brief Returns the number of tetrahedra of the mesh
   * \return numTets the number of tetrahedra.
   */
  int getNumTets() const { return m_num_tets; };

  /*!
   * \brief Clears all internal data-structures
   */
  void clear();

  /*!
   * \brief Reads in the mesh from a Pro/E file
   * \pre path to input file has been set by calling `setFileName()`
   * \return status set to zero on success; set to a non-zero value otherwise.
   */
  virtual int read();

  /*!
   * \brief Sets the bounding box to include nodes.
   * 
   * When reading in ProE tets, include only those tets that are fully
   * contained in the bounding box.  If the box is left unset, or if set
   * to an invalid box, all tets are included.
   */
  virtual void setBoundingBox(BBox3D& box);

  /*!
   * \brief Sets the bounding box to include nodes.
   * 
   * When reading in ProE tets, include only those tets that have at least 
   * one node in the bounding box.  If the box is left unset, or if set
   * to an invalid box, all tets are included.
   */
  virtual void setInclusiveBoundingBox(BBox3D& box);

  /*!
   * \brief Specifies a test to retain a tetrahedron.
   * 
   * When reading in ProE tets, include only those tets for which p returns
   * true.
   */
  virtual void setTetPred(const TetPred& p) { m_tetPredicate = p; }

  /*!
   * \brief Stores the Pro/E data in the supplied unstructured mesh object.
   * \param [in,out] mesh pointer to the unstructured mesh.
   * \pre mesh != nullptr.
   */
  void getMesh(mint::UnstructuredMesh<mint::SINGLE_SHAPE>* mesh);

protected:
  std::string m_fileName;

  axom::IndexType m_num_nodes;
  axom::IndexType m_num_tets;

  std::vector<double> m_nodes;
  std::vector<int> m_tets;

  TetPred m_tetPredicate;

  /*!
   * \brief Compact internal mesh storage arrays
   * \param [in] retain_vertex Flag of all vertices to retain after compaction
   * \param [in] elt_count Number of elements stored in array
   * 
   * This method compacts the vertex array and resizes the vertex and the
   * element arrays to their minimum sizes.
   */
  void compact_arrays(slam::BitSet& retain_vertex, int elt_count);

private:
  DISABLE_COPY_AND_ASSIGNMENT(ProEReader);
  DISABLE_MOVE_AND_ASSIGNMENT(ProEReader);
};

}  // namespace quest
}  // namespace axom

#endif  // QUEST_PROEREADER_HPP_
