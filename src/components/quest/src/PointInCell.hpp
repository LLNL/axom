/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef QUEST_POINT_IN_CELL_HPP_
#define QUEST_POINT_IN_CELL_HPP_


#include "axom/config.hpp"
#include "axom/Macros.hpp"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#include "quest/ImplicitGrid.hpp"

#include "mint/UnstructuredMesh.hpp"
#include "mint/vtk_utils.hpp"
#include "mint/FieldData.hpp"
#include "mint/FieldVariable.hpp"


#ifdef AXOM_USE_MFEM
# include "quest/PointInCell_impl_mfem.hpp"
#endif

#ifndef MFEM_ISO_TRANS_HAS_DEBUG_STRUCT
//#  define MFEM_ISO_TRANS_HAS_DEBUG_STRUCT
#endif


namespace axom {
namespace quest {

struct quest_point_in_cell_mint_tag {};

namespace detail {

/*!
 * \class
 * \brief A traits class for the mesh associated with a PointInCell query.
 *
 * An implementation of PointInCellTraits must define:
 * \arg A MeshType (e.g. mfem::Mesh or axom::mint::Mesh)
 * \arg An IndexType (e.g. int)
 * \arg An IndexType variable named NO_CELL (e.g. with value -1), indicating an
 *   invalid index in the mesh
 */
template<typename quest_point_in_cell_mesh_tag>
struct PointInCellTraits;

/*!
 * \class
 * \brief A wrapper for access to the mesh instance associated with a PointInCell query
 *
 * An implementation of a PointInCellMeshWrapper must expose the following API:
 * \arg computeBoundingBox(...)
 * \arg findInSpace(...)
 * \arg getIsoparametricCoords(...)
 * \arg numElements()
 * \arg meshDimension()
 *
 * \TODO -- fill this in and ensure it is accurate
 */
template <typename quest_point_in_cell_mesh_tag>
class PointInCellMeshWrapper;

template<>
struct PointInCellTraits<quest_point_in_cell_mint_tag>
{
  typedef axom::mint::Mesh MeshType;
  typedef int IndexType;

  static const IndexType NO_CELL;
};

/*!
 * \class
 *
 * \brief A class to encapsulate locating a point within a computational mesh
 */
template<int NDIMS, typename mesh_tag>
class PointFinder
{
public:
  typedef quest::ImplicitGrid<NDIMS> GridType;

  typedef typename GridType::SpacePoint SpacePoint;
  typedef typename GridType::SpatialBoundingBox SpatialBoundingBox;

  typedef PointInCellMeshWrapper<mesh_tag> MeshWrapperType;
  typedef typename MeshWrapperType::IndexType IndexType;

public:
  PointFinder(MeshWrapperType* meshWrapper, int* res)
    : m_meshWrapper(meshWrapper)
  {
    SLIC_ASSERT( m_meshWrapper != AXOM_NULLPTR);

    const int numElem = m_meshWrapper->numElements();

    // setup bounding boxes -- Slightly scaled for robustness

    SpatialBoundingBox meshBBox;
    m_eltBBoxes = std::vector<SpatialBoundingBox>(numElem);
    double EPS = 1e-8;
    double bboxScaleFactor = 1 + EPS;
    m_meshWrapper->template computeBoundingBoxes<NDIMS>(bboxScaleFactor, m_eltBBoxes, meshBBox);
    //SLIC_INFO("Mesh bounding box is: " << meshBBox );

    // initialize implicit grid
    typedef axom::primal::Point<int, NDIMS> GridResolution;
    GridResolution gridRes;
    if(res != AXOM_NULLPTR)
    {
      gridRes = GridResolution(res);
    }
    m_grid.initialize(meshBBox, gridRes, numElem);

    // add mesh elements to grid
    for(int i=0; i< numElem; ++i)
    {
      m_grid.insert( m_eltBBoxes[i], i);
    }
  }

  IndexType locatePoint(const double* pos, double* isoparametric)
  {
    typedef typename GridType::BitsetType BitsetType;

    SLIC_ASSERT( pos != AXOM_NULLPTR);
    SpacePoint pt(pos);


    IndexType containingCell = PointInCellTraits<mesh_tag>::NO_CELL;

    SpacePoint isopar;

    // Note: ImplicitGrid::getCandidates() checks the mesh bounding box
    BitsetType candidates = m_grid.getCandidates(pt);

    // SLIC_INFO("Candidates for space point " <<  pt);
    bool foundContainingCell = false;
    for(std::size_t eltIdx = candidates.find_first();
        !foundContainingCell && eltIdx != BitsetType::npos;
        eltIdx = candidates.find_next( eltIdx) )
    {
      // First check that pt is in bounding box of element
      if( m_eltBBoxes[eltIdx].contains(pt) )
      {
        // if isopar is in the proper range
        if( m_meshWrapper->getIsoparametricCoords(eltIdx, pt.data(), isopar.data() ) )
        {
          // then we have found the cellID
          foundContainingCell = true;
          containingCell = eltIdx;
        }
      }
    }

    // Copy data back to input parameter isoparametric, if necessary
    if(isoparametric != AXOM_NULLPTR)
    {
      isopar.array().to_array(isoparametric);
    }

    return containingCell;
  }

//  void enableDebugMeshGeneration()
//  {
//    if(m_queryPathsMeshDumper == AXOM_NULLPTR)
//    {
//      m_queryPathsMeshDumper = new detail::PICQueryMeshDumper<NDIMS>;
//    }
//  }
//  void printDebugMesh(const std::string& filename)
//  {
//    if( m_queryPathsMeshDumper != AXOM_NULLPTR)
//    {
//      m_queryPathsMeshDumper->outputDebugMesh(filename);
//    }
//  }
//  detail::PICQueryMeshDumper<NDIMS>* m_queryPointsMeshDumper;
//  detail::PICQueryMeshDumper<NDIMS>* m_queryPathsMeshDumper;

private:
  GridType    m_grid;

  MeshWrapperType* m_meshWrapper;

  std::vector<SpatialBoundingBox> m_eltBBoxes;
};

/*!
 * Utility class for visualizing the point in cell queries.
 *
 * \note This class can output a more detailed mesh if mfem's Isoparametric
 * Transformation class is annotated with the following internal struct:
 * \code
 *    struct TransformBackDebug {
 *       int numIters;
 *       DenseMatrix points;
 *    } tbDebug;
 *  \endcode
 * This is not included in the mfem release:
 */
template<int NDIMS>
class PICQueryMeshDumper
{
public:
  typedef quest::ImplicitGrid<NDIMS> GridType;
  typedef typename GridType::SpacePoint SpacePoint;

public:
  PICQueryMeshDumper() : m_transformBackDebugMesh(NDIMS), m_numQueries(0) {}

public:

  /*!
   * \brief Adds the path of the query for point \a pt associated with
   * transformation \a tr to the mesh
   *
   * \param [in] pt The query point
   * \param [in] isopar The isoparametric coordianates found by the query
   * \param [in] status The status of the query
   * \param [in] tr The tranformation from the query.
   * \note Parameter \a tr is an [in] parameter, but cannot be \a const
   * since this function calls a non-const member function of \a tr
   */
  void addQueryPoint(const SpacePoint& pt, const SpacePoint& isopar, int status, mfem::IsoparametricTransformation& tr)
  {

    int meshStartIndex = m_transformBackDebugMesh.getMeshNumberOfNodes();

    // Transform last found reference pt into space
    mfem::Vector y;
    mfem::IntegrationPoint ipRef;
    ipRef.Set(isopar.data(), NDIMS);
    tr.Transform(ipRef, y);

  #ifdef MFEM_ISO_TRANS_HAS_DEBUG_STRUCT
    //SLIC_INFO("Used " << tr.tbDebug.numIters << " iterations." );

    tr.tbDebug.points.SetCol( tr.tbDebug.numIters, y );

    // Insert the segment vertex positions into the mesh
    for(int i=0; i< tr.tbDebug.numIters; ++i)
    {
      double ptx = tr.tbDebug.points(0, i);
      double pty = tr.tbDebug.points(1, i);
      if(NDIMS==2)
      {
        m_transformBackDebugMesh.insertNode(ptx,pty);
      }
      else if(NDIMS == 3)
      {
        double ptz = tr.tbDebug.points(2, i);
        m_transformBackDebugMesh.insertNode(ptx,pty,ptz);
      }
    }

    // Add segments to mesh and record associated 'field' data
    for(int i=0; i< tr.tbDebug.numIters-1; ++i)
    {
      int currVertIdx = meshStartIndex + i;
      int indices[2] = { currVertIdx, currVertIdx+1 };
      m_transformBackDebugMesh.insertCell(indices, MINT_SEGMENT, 2);

      m_query_spacept.push_back( pt);
      m_query_refpt.push_back(isopar);
      m_query_index.push_back( m_numQueries );
      m_query_return_code.push_back( status );
      m_query_num_iters.push_back( tr.tbDebug.numIters );
      m_query_iter.push_back(i);
    }
#else

    // Add the point the mesh
    m_transformBackDebugMesh.insertNode( y.GetData() );

    // Add the fields: status, pt, isopt, index
    int indices[1] = { meshStartIndex };
    m_transformBackDebugMesh.insertCell( indices, MINT_VERTEX, 1);

    m_query_spacept.push_back( pt);
    m_query_refpt.push_back(isopar);
    m_query_index.push_back( m_numQueries );
    m_query_return_code.push_back( status );

#endif

    m_numQueries++;
  }

  /*! Dumps the mesh to disk */
  void outputDebugMesh(const std::string& filename)
  {
    int size = m_query_index.size();

    if(size == 0)
    {
        SLIC_INFO("No mesh data.");
        return;
    }

    // Add integer fields
    int* fIndex = addIntField("query_index",size);
    int* fStatus = addIntField("query_status",size);

    bool hasNIters = m_query_num_iters.size() > 0;
    int* fNIters = hasNIters ?addIntField("query_num_iters",size) : AXOM_NULLPTR;

    bool hasIters = m_query_iter.size() > 0;
    int* fIter = hasIters ? addIntField("query_iter",size) : AXOM_NULLPTR;

    for(int i=0; i < size; ++i)
    {
      fIndex[i]  = m_query_index[i];
      fStatus[i] = m_query_return_code[i];

      if(hasNIters)
        fNIters[i] = m_query_num_iters[i];
      if(hasIters)
        fIter[i]   = m_query_iter[i];
    }

    // Add (vector-like) point fields
    double* fPt[3];
    fPt[0] = addDoubleField("query_pt/x", size);
    fPt[1] = addDoubleField("query_pt/y", size);
    fPt[2] = NDIMS ==3 ? addDoubleField("query_pt/z", size) : AXOM_NULLPTR;

    double* fIsoPt[3];
    fIsoPt[0] = addDoubleField("query_iso_pt/x", size);
    fIsoPt[1] = addDoubleField("query_iso_pt/y", size);
    fIsoPt[2] = NDIMS ==3 ? addDoubleField("query_iso_pt/z", size) : AXOM_NULLPTR;

    for(int i=0; i < size; ++i)
    {
      for(int d=0; d < NDIMS; ++d)
      {
        fPt[d][i] = m_query_spacept[i][d];
        fIsoPt[d][i] = m_query_refpt[i][d];
      }
    }

    // Output mesh to disk
    mint::write_vtk(&m_transformBackDebugMesh, filename);
  }

private:

  /*! Utility function to add an integer field to the mint mesh */
  int* addIntField(const std::string& name, int size)
  {
    mint::FieldData* CD = m_transformBackDebugMesh.getCellFieldData();
    CD->addField( new mint::FieldVariable< int >(name, size) );

    int* fld = CD->getField( name )->getIntPtr();

    SLIC_ASSERT(fld != AXOM_NULLPTR);
    return fld;
  }

  /*! Utility function to add a double field to the mint mesh */
  double* addDoubleField(const std::string& name, int size)
  {
    mint::FieldData* CD = m_transformBackDebugMesh.getCellFieldData();
    CD->addField( new mint::FieldVariable< double >(name, size) );

    double* fld = CD->getField( name )->getDoublePtr();

    SLIC_ASSERT(fld != AXOM_NULLPTR);
    return fld;
  }

private:

  mint::UnstructuredMesh<MINT_MIXED_CELL> m_transformBackDebugMesh;

  std::vector<int> m_query_index;
  std::vector<int> m_query_return_code;
  std::vector<int> m_query_num_iters;
  std::vector<int> m_query_iter;
  std::vector<SpacePoint> m_query_spacept;
  std::vector<SpacePoint> m_query_refpt;

  int m_numQueries;
};


} // end namespace detail

/*!
 * \class PointInCell
 *
 * \brief A spatial querying class to accelerate Point-In-Cell queries.
 *
 * A point in cell query over a given mesh determines the index of the mesh cell
 * containing a given element.  Additionally, when an element is found, the query
 * can find the isoparametric coordinates of the point within the element, and the
 * interpolation weights of the point w.r.t. the element's basis functions.
 *
 * \note This class was designed to support point in cell queries against
 * a 2D or 3D computational mesh. The queries to the mesh are wrapped in a PointInCellMeshWrapper
 * class templated on a mesh tag. To extend this design, one must create a new mesh tag, e.g.
 * and add custom implementations of PointInCellMeshWrapper<custom_mesh_tag>
 * and PointInCellTraits<custom_mesh_tag> in the axom::quest::detail namespace.
 * See PointInCell_impl_mfem.hpp for an example.
 */
template<typename mesh_tag>
class PointInCell
{
public:
  typedef detail::PointInCellMeshWrapper<mesh_tag> MeshWrapperType;
  typedef typename detail::PointInCellTraits<mesh_tag>::MeshType MeshType;
  typedef int IndexType;

  typedef detail::PointFinder<2, mesh_tag> PointFinder2D;
  typedef detail::PointFinder<3, mesh_tag> PointFinder3D;


  /*!
   * Construct a point in cell query structure over a given mfem mesh
   * \param mesh A pointer to an mfem mesh
   * \param resolution If the resolution is not provided, we use a heuristic to set the resolution
   */
  PointInCell(MeshType* mesh, int* resolution = AXOM_NULLPTR)
    : m_meshWrapper(mesh), m_pointFinder2D(AXOM_NULLPTR), m_pointFinder3D(AXOM_NULLPTR)
  {
    // Allocate a 2D or 3D PointFinder, depending on the mesh dimension.
    // Note: Only one of these will be allocated in a PointInCell instance

    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      m_pointFinder2D = new PointFinder2D(&m_meshWrapper, resolution);
      break;
    case 3:
      m_pointFinder3D = new PointFinder3D(&m_meshWrapper, resolution);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }
  }

  ~PointInCell()
  {
    if( m_pointFinder2D != AXOM_NULLPTR)
    {
      delete m_pointFinder2D;
      m_pointFinder2D = AXOM_NULLPTR;
    }

    if( m_pointFinder3D != AXOM_NULLPTR)
    {
      delete m_pointFinder3D;
      m_pointFinder3D = AXOM_NULLPTR;
    }
  }

  /*! Attempt to find the index of the cell containing the given point.
   *  If found, also return isoparametric coords in out param isopar
   */
  /*!
   * \brief Attempts to find the index of the mesh cell containing query point pt
   * \param[in] pt The query point
   * \param[out] isoparametric When pt is within a cell, this point returns
   *             the isoparametric coordinates of pt within the cell.
   * \return The index of the mesh cell containing pt, when one exists,
   *         otherwise returns the special value PointInCell::NO_CELL
   */
   
  IndexType locatePoint(const double* pos, double* isopar)
  {
    IndexType cellIndex = detail::PointInCellTraits<mesh_tag>::NO_CELL;

    switch(m_meshWrapper.meshDimension())
    {
    case 2:
      cellIndex = m_pointFinder2D->locatePoint(pos, isopar);
      break;
    case 3:
      cellIndex = m_pointFinder3D->locatePoint(pos, isopar);
      break;
    default:
      SLIC_ERROR("Point in Cell query only defined for 2D or 3D meshes.");
      break;
    }

    return cellIndex;
  }

  /** Attempt to find the isoparametric coords of a point pos in the given mesh cell */
  bool locatePointInCell(IndexType cellIdx, const double* pos, double* isopar)
  {
    // TODO: Check against element's bbox

    return m_meshWrapper.getIsoparametricCoords(cellIdx, pos, isopar);
  }

  /*! Transform from isoparametric space to physical space */
  /*!
   * Returns the point in space corresponding to isoparameteric coordinates isopar
   * of the mesh element with index eltIdx
   */
  
  void findInSpace(IndexType cellIdx, const double* isopar, double* pos)
  {
    m_meshWrapper.findInSpace(cellIdx, isopar, pos);
  }

private:
  MeshWrapperType m_meshWrapper;

  PointFinder2D* m_pointFinder2D;
  PointFinder3D* m_pointFinder3D;
  
  // double m_eps; // for scaling the bounding boxes? 
};




} // end namespace quest
} // end namespace axom

#endif // QUEST_POINT_IN_CELL_HPP_
