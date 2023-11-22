// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file MarchingCubes.hpp
 *
 * \brief Consists of classes implementing marching cubes algorithm to
 * compute isocontour from a scalar field in a blueprint mesh.
 */

#ifndef AXOM_QUEST_MARCHINGCUBES_H_
#define AXOM_QUEST_MARCHINGCUBES_H_

#include "axom/config.hpp"

// Implementation requires Conduit.
#ifdef AXOM_USE_CONDUIT

  // Axom includes
  #include "axom/mint/mesh/UnstructuredMesh.hpp"

  // Conduit includes
  #include "conduit_node.hpp"

  // C++ includes
  #include <string>

  // Add some helper preprocessor defines for using OPENMP, CUDA, and HIP policies
  // within the marching cubes implementation.
  #if defined(AXOM_USE_RAJA)
    #ifdef AXOM_USE_OPENMP
      #define _AXOM_MARCHINGCUBES_USE_OPENMP
    #endif
    #if defined(AXOM_USE_CUDA) && defined(AXOM_USE_UMPIRE)
      #define _AXOM_MARCHINGCUBES_USE_CUDA
    #endif
    #if defined(AXOM_USE_HIP) && defined(AXOM_USE_UMPIRE)
      #define _AXOM_MARCHINGCUBES_USE_HIP
    #endif
  #endif

namespace axom
{
namespace quest
{
/*!
  @brief Enum for runtime execution policy

  The policy implicitly selects the execution space and allocator id.
*/
enum class MarchingCubesRuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2,
  hip = 3
};

/// Utility function to allow formating a MarchingCubesRuntimePolicy
inline auto format_as(MarchingCubesRuntimePolicy pol)
{
  return fmt::underlying(pol);
}

namespace detail
{
namespace marching_cubes
{
template <int DIM, typename ExecSpace, typename SequentialLoopPolicy>
class MarchingCubesImpl;
}
}  // namespace detail

class MarchingCubesSingleDomain;

/*!
 * \@brief Class implementing marching cubes to compute a contour
 * mesh from a scalar function on an input mesh.
 *
 * This implementation is for the original 1987 algorithm:
 * Lorensen, William E.; Cline, Harvey E. (1 August 1987).
 * "Marching cubes: A high resolution 3D surface construction algorithm".
 * ACM SIGGRAPH Computer Graphics. 21 (4): 163-169
 *
 * Implementation is for 2D (marching squares) and 3D (marching
 * cubes).
 *
 * The input mesh is a Conduit::Node following the Mesh Blueprint
 * convention.  The mesh must be in multi-domain format.  For
 * single-domain, see MarchingCubesSingleDomain.
 *
 * Usage example:
 * @beginverbatim
 *   void foo( conduit::Node &meshNode,
 *             const std::string &coordsName,
 *             const std::string &functionName,
 *             double contourValue )
 *   {
 *     MarchingCubes mc(meshNode, coordsName);
 *     setFunctionField(functionName);
 *     mc.computeIsocontour(contourValue);
 *     axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE>
 *       contourMesh(3, min::CellType::Triangle);
 *     mc.populateContourMesh( contourMesh, "cellIdField");
 *   }
 * @endverbatim
 *
 * To avoid confusion between the two meshes, we refer to the input
 * mesh with the scalar function as "parent" and the generated mesh
 * as the "contour".
 */
class MarchingCubes
{
public:
  using RuntimePolicy = MarchingCubesRuntimePolicy;
  /*!
   * \brief Constructor sets up computational mesh and data for running the
   * marching cubes algorithm.
   *
   * \param [in] runtimePolicy A value from MarchingCubesRuntimePolicy.
   *             The simplest policy is RuntimePolicy::seq, which specifies
   *             running sequentially on the CPU.
   * \param [in] bpMesh Blueprint multi-domain mesh containing scalar field.
   * \param [in] topologyName Name of Blueprint topology to use in \a bpMesh.
   * \param [in] maskField Cell-based std::int32_t mask field.  If provided,
   *             cells where this field evaluates to false are skipped.
   *
   * Array data in \a dom must be accessible in the \a runtimePolicy
   * environment.  It's an error if not, e.g., using CPU memory with
   * a GPU policy.
   *
   * Some data from \a bpMesh may be cached by the constructor.
   * Any change to it after the constructor leads to undefined behavior.
   *
   * The mesh coordinates should be contiguous.  See
   * conduit::blueprint::is_contiguous().  In the future, this
   * requirement may be relaxed, possibly at the cost of a
   * transformation and storage of the temporary contiguous layout.
   */
  MarchingCubes(RuntimePolicy runtimePolicy,
                const conduit::Node &bpMesh,
                const std::string &topologyName,
                const std::string &maskField = {});

  /*!
    @brief Set the field containing the nodal function.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void setFunctionField(const std::string &fcnField);

  /*!
   \brief Computes the isocontour.
   \param [in] contourVal isocontour value
   */
  void computeIsocontour(double contourVal = 0.0);

  //!@brief Get number of cells in the generated contour mesh.
  axom::IndexType getContourCellCount() const;

  //!@brief Get number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const;

  /*!
    @brief Put generated contour in a mint::UnstructuredMesh.
    @param mesh Output contour mesh
    @param cellIdField Name of field to store the (axom::IndexType)
      parent cell ids. If omitted, the data is not provided.
    @param domainIdField Name of field to store the (axom::IndexType)
      parent domain ids. If omitted, the data is not provided.

    If the fields aren't in the mesh, they will be created.

    Blueprint allows users to specify ids for the domains.  If
    "state/domain_id" exists in the domains, it is used as the domain
    id.  Otherwise, the domain's interation index within the
    multidomain mesh is used.
  */
  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
    const std::string &cellIdField = {},
    const std::string &domainIdField = {});

private:
  MarchingCubesRuntimePolicy m_runtimePolicy;

  //! @brief Single-domain implementations.
  axom::Array<std::unique_ptr<MarchingCubesSingleDomain>> m_singles;
  const std::string m_topologyName;
  std::string m_fcnFieldName;
  std::string m_fcnPath;
  std::string m_maskFieldName;
  std::string m_maskPath;

  void setMesh(const conduit::Node &bpMesh);
};

/*!
 * \@brief Class implementing marching cubes algorithm for a single
 *  domain.
 *
 * \sa MarchingCubes
 */
class MarchingCubesSingleDomain
{
  template <int DIM, typename ExecSpace, typename SequentialLoopPolicy>
  friend class detail::marching_cubes::MarchingCubesImpl;

public:
  using RuntimePolicy = MarchingCubesRuntimePolicy;
  /*!
   * \brief Constructor for applying algorithm in a single domain.
   * See MarchingCubes for the multi-domain implementation.
   *
   * \param [in] runtimePolicy A value from MarchingCubesRuntimePolicy.
   *             The simplest policy is RuntimePolicy::seq, which specifies
   *             running sequentially on the CPU.
   * \param [in] dom Blueprint single-domain mesh containing scalar field.
   * \param [in] topologyName Name of Blueprint topology to use in \a dom
   * \param [in] maskField Cell-based std::int32_t mask field.  If provided,
   *             cells where this field evaluates to false are skipped.
   *
   * Array data in \a dom must be accessible in the \a runtimePolicy
   * environment.  It's an error if not, e.g., using CPU memory with
   * a GPU policy.
   *
   * Some data from \a dom may be cached by the constructor.
   * Any change to it after the constructor leads to undefined behavior.
   * See setDomain(const conduit::Node &) for requirements on \a dom.
   *
   * The mesh coordinates should be stored contiguously.  See
   * conduit::blueprint::is_contiguous().  In the future, this
   * requirement may be relaxed, possibly at the cost of a
   * transformation and storage of the temporary contiguous layout.
   */
  MarchingCubesSingleDomain(RuntimePolicy runtimePolicy,
                            const conduit::Node &dom,
                            const std::string &topologyName,
                            const std::string &maskfield);

  /*!
    @brief Specify the field containing the nodal scalar function
    in the input mesh.
    \param [in] fcnField Name of node-based scalar function values.
  */
  void setFunctionField(const std::string &fcnField);

  /*!
    @brief Get the Blueprint domain id specified in \a state/domain_id
    if it is provided, or use the given default if not provided.
  */
  int getDomainId(int defaultId) const;

  /*!
   * \brief Compute the isocontour.
   *
   * \param [in] contourVal isocontour value
   *
   * Compute isocontour using the marching cubes algorithm.
   * To get the isocontour mesh, use populateContourMesh.
   */
  void computeIsocontour(double contourVal = 0.0);

  //!@brief Get number of cells in the generated contour mesh.
  axom::IndexType getContourCellCount() const
  {
    SLIC_ASSERT_MSG(
      m_impl,
      "There is no contour mesh until you call computeIsocontour()");
    axom::IndexType cellCount = m_impl->getContourCellCount();
    return cellCount;
  }

  //!@brief Get number of nodes in the generated contour mesh.
  axom::IndexType getContourNodeCount() const
  {
    return m_ndim * getContourCellCount();
  }

  /*!
    @brief Put generated contour surface in a mint::UnstructuredMesh
    object.

    @param mesh Output mesh
    @param cellIdField Name of field to store the prent cell ids.
      If omitted, the data is not copied.
  */
  void populateContourMesh(
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
    const std::string &cellIdField = {}) const
  {
    m_impl->populateContourMesh(mesh, cellIdField);
  }

  /*!
    @brief Determine whether a given \a MarchingCubesRuntimePolicy is
    valid for the Axom build configuration.
  */
  inline bool isValidRuntimePolicy(MarchingCubesRuntimePolicy policy) const
  {
    switch(policy)
    {
    case MarchingCubesRuntimePolicy::seq:
      return true;

    case MarchingCubesRuntimePolicy::omp:
  #ifdef _AXOM_MARCHINGCUBES_USE_OPENMP
      return true;
  #else
      return false;
  #endif

    case MarchingCubesRuntimePolicy::cuda:
  #ifdef _AXOM_MARCHINGCUBES_USE_CUDA
      return true;
  #else
      return false;
  #endif
    case MarchingCubesRuntimePolicy::hip:
  #ifdef _AXOM_MARCHINGCUBES_USE_HIP
      return true;
  #else
      return false;
  #endif
    }

    return false;
  }

private:
  MarchingCubesRuntimePolicy m_runtimePolicy;
  /*!
    \brief Computational mesh as a conduit::Node.
  */
  const conduit::Node *m_dom;
  int m_ndim;

  //!@brief Name of Blueprint topology in m_dom.
  const std::string m_topologyName;

  std::string m_fcnFieldName;
  //!@brief Path to nodal scalar function in m_dom.
  std::string m_fcnPath;

  const std::string m_maskFieldName;
  //!@brief Path to mask in m_dom.
  const std::string m_maskPath;

  /*!
    @brief Base class for implementations templated on dimension
    and execution space.

    This class allows m_impl to refer to any implementation used
    at runtime.
  */
  struct ImplBase
  {
    //!@brief Prepare internal data for operating on the given domain.
    virtual void initialize(const conduit::Node &dom,
                            const std::string &topologyName,
                            const std::string &fcnPath,
                            const std::string &maskPath) = 0;
    //!@brief Set the contour value
    virtual void setContourValue(double contourVal) = 0;
    //@{
    //!@name Phases of the computation
    //!@brief Mark domain cells that cross the contour.
    virtual void markCrossings() = 0;
    //!@brief Precompute some metadata for contour mesh.
    virtual void scanCrossings() = 0;
    //!@brief Generate the contour mesh in internal data format.
    virtual void computeContour() = 0;
    //!@brief Get the number of contour mesh cells generated.
    //@}
    virtual axom::IndexType getContourCellCount() const = 0;
    /*!
      @brief Populate output mesh object with generated contour.

      Note: Output format is in flux.  We will likely output
      a blueprint object in the future.
    */
    virtual void populateContourMesh(
      axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> &mesh,
      const std::string &cellIdField) const = 0;
    virtual ~ImplBase() { }
  };

  std::unique_ptr<ImplBase> m_impl;
  //!@brief Allocate implementation object and set m_impl.
  void allocateImpl();

  /*!
   * \brief Set the blueprint single-domain mesh.
   *
   * Some data from \a dom may be cached by the constructor.
   */
  void setDomain(const conduit::Node &dom);

};  // class MarchingCubesSingleDomain

}  // namespace quest
}  // namespace axom

#endif  // AXOM_USE_CONDUIT
#endif  // AXOM_QUEST_MARCHINGCUBES_H_
