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


#ifndef AXOM_USE_MFEM
#  error Mfem is a required dependency for quest::PointInCell.
#endif

#ifndef MFEM_ISO_TRANS_HAS_DEBUG_STRUCT
//#  define MFEM_ISO_TRANS_HAS_DEBUG_STRUCT
#endif

#include "mfem.hpp"

namespace axom {
namespace quest {

namespace detail {

template < int NDIMS, typename MeshType > class PICMeshWrapper;

/*!
 * \brief Wraps MFEM mesh access functionality for the PointInCell class
 */
template < int NDIMS>
class PICMeshWrapper<NDIMS, mfem::Mesh>
{
public:
  typedef axom::primal::Point<double, NDIMS> SpacePoint;
  typedef axom::primal::BoundingBox<double, NDIMS> SpatialBoundingBox;
  typedef int IndexType;

  PICMeshWrapper(mfem::Mesh* mesh) : m_mesh(mesh), m_isInitialized(false)
  {
    // Some sanity checks
    SLIC_ASSERT( m_mesh != AXOM_NULLPTR);
    SLIC_ASSERT( m_mesh->Dimension() == NDIMS);

    m_isHighOrder = (m_mesh->GetNodalFESpace() != AXOM_NULLPTR) && (m_mesh->GetNE() > 0);

    // Initialize bounding boxes -- split into two functions for simplicity
    double const EPS = 1e-8;
    double bboxScaleFactor = 1. + EPS;
    const int numMeshElements = numElements();
    m_elementBoundingBox = std::vector<SpatialBoundingBox>(numMeshElements);

    if( m_isHighOrder)
    {
      initializeBoundingBoxesHighOrder(bboxScaleFactor);
    }
    else
    {
      initializeBoundingBoxesLowOrder(bboxScaleFactor);
    }

    m_isInitialized = true;
  }

  /*! Predicate to check if the given basis is positive (i.e. Bernstein) */
  static bool isPositiveBasis(const mfem::FiniteElementCollection* fec)
  {
    // HACK -- is this sufficient?  Is there a better way to do this?

    if(fec == AXOM_NULLPTR)
    {
      return false;
    }

    if(const mfem::H1_FECollection* h1Fec = dynamic_cast<const mfem::H1_FECollection*>(fec))
    {
      return h1Fec->GetBasisType() == mfem::BasisType::Positive;
    }

    if(const mfem::L2_FECollection* l2Fec = dynamic_cast<const mfem::L2_FECollection*>(fec))
    {
      return l2Fec->GetBasisType() == mfem::BasisType::Positive;
    }

    if( dynamic_cast<const mfem::NURBSFECollection*>(fec)       ||
        dynamic_cast<const mfem::LinearFECollection*>(fec)      ||
        dynamic_cast<const mfem::QuadraticPosFECollection*>(fec) )
    {
      return true;
    }

    return false;
  }

  /*!
   * \brief Utility function to get a positive (i.e. Bernstein) collection of bases
   * corresponding to the given FiniteElementCollection.
   *
   * \return A pointer to a newly allocated FiniteElementCollection
   * corresponding to \a fec
   * \note   It is the user's responsibility to deallocate this pointer.
   * \pre    \a fec is not already positive
   */
  static mfem::FiniteElementCollection* getCorrespondingPositiveFEC(const mfem::FiniteElementCollection* fec, int order, int mapType)
  {
    //   NOTE(KW) Should this use a map<string, string> instead ?
    //         It could then use FECollection::New(string)

   // If fec is already positive, return it?
    SLIC_CHECK_MSG( !isPositiveBasis(fec),
        "This function is only meant to be called on non-positive finite element collection"
        );

    // Attempt to find the corresponding positive H1 fec
    if(const mfem::H1_FECollection* h1Fec = dynamic_cast<const mfem::H1_FECollection*>(fec))
    {
      return new mfem::H1_FECollection(order, NDIMS, mfem::BasisType::Positive);
    }

    // Attempt to find the corresponding positive L2 fec
    if(const mfem::L2_FECollection* l2Fec = dynamic_cast<const mfem::L2_FECollection*>(fec))
    {
      // should we throw a not supported error here?
      return new mfem::L2_FECollection(order, NDIMS, mfem::BasisType::Positive, mapType);
    }

    // Attempt to find the corresponding quadratic or cubic fec
    if(dynamic_cast<const mfem::QuadraticFECollection*>(fec) ||
       dynamic_cast<const mfem::CubicFECollection*>(fec) )
    {
      SLIC_ASSERT( order == 2 || order == 3);
      return new mfem::H1_FECollection(order, NDIMS, mfem::BasisType::Positive);
    }

    // Give up -- return NULL
    return AXOM_NULLPTR;
  }


  /*! Get the number of elements in the mesh */
  int numElements() const { return m_mesh->GetNE(); }

  /*! Get a pointer to the mesh */
  mfem::Mesh* getMesh() const { return m_mesh; }

  /*! Return the bounding box of the mesh */
  const SpatialBoundingBox& meshBoundingBox() const
  {
    return m_meshBoundingBox;
  }

  /*! Return the bounding box of mesh element with index idx */
  const SpatialBoundingBox& elementBoundingBox(IndexType idx) const
  {
    SLIC_ASSERT( 0 <= idx && idx < numElements() );
    return m_elementBoundingBox[idx];
  }

private:
  /*!
   * Helper function to initialize the bounding boxes for a high order mfem mesh
   *
   * \param bboxScaleFactor Scale factor for expanding bounding boxes
   * \note This function can only be called before the class is initialized
   * \pre This function can only be called when m_isHighOrder is true
   * \pre bboxScaleFactor must be greater than or equal to 1.
   */
  void initializeBoundingBoxesHighOrder(double bboxScaleFactor)
  {
    // Sanity checks
    SLIC_ASSERT( !m_isInitialized );
    SLIC_ASSERT( m_isHighOrder );    SLIC_ASSERT( !m_isInitialized );
    SLIC_ASSERT(bboxScaleFactor >= 1.);

    /// Generate (or access existing) positive (Bernstein) nodal grid function
    const mfem::FiniteElementSpace* nodalFESpace = m_mesh->GetNodalFESpace();

    mfem::GridFunction* positiveNodes = AXOM_NULLPTR;
    bool mustDeleteNodes = false;

    const mfem::FiniteElementCollection* nodalFEColl = nodalFESpace->FEColl();

    // Check if grid function is positive, if not create positive grid function
    if( PICMeshWrapper::isPositiveBasis( nodalFEColl ) )
    {
      positiveNodes = m_mesh->GetNodes();
    }
    else
    {
      // Assume that all elements of the mesh have the same order and geom type
      int order = nodalFESpace->GetOrder(0);
      int GeomType = m_mesh->GetElementBaseGeometry(0);
      int mapType = (nodalFEColl != AXOM_NULLPTR)
          ? nodalFEColl->FiniteElementForGeometry(GeomType)->GetMapType()
          : static_cast<int>(mfem::FiniteElement::VALUE);

      mfem::FiniteElementCollection* posFEColl =
          PICMeshWrapper::getCorrespondingPositiveFEC(nodalFEColl, order, mapType);

      SLIC_ASSERT_MSG(posFEColl != AXOM_NULLPTR,
          "Problem generating a positive finite element collection "
          << "corresponding to the mesh's '"<< nodalFEColl->Name()
          << "' finite element collection.");

      if(posFEColl != AXOM_NULLPTR)
      {
        // Create a positive (Bernstein) grid function for the nodes
        mfem::FiniteElementSpace* posFESpace =
            new mfem::FiniteElementSpace(m_mesh, posFEColl, NDIMS);
        positiveNodes = new mfem::GridFunction(posFESpace);

        // m_bernsteinNodes takes ownership of posFEColl's memory
        positiveNodes->MakeOwner(posFEColl);
        mustDeleteNodes = true;

        // Project the nodal grid function onto this
        positiveNodes->ProjectGridFunction(*(m_mesh->GetNodes()));
      }
    }

    // Output some information
    SLIC_INFO("Mesh nodes fec -- " << nodalFEColl->Name()
        << " with ordering " << nodalFESpace->GetOrdering()
        << "\n\t -- Positive nodes are fec -- "
        << positiveNodes->FESpace()->FEColl()->Name()
        << " with ordering " << positiveNodes->FESpace()->GetOrdering() );


    /// For each element, compute bounding box, and overall mesh bbox
    mfem::Array<int> dofIndices;
    mfem::FiniteElementSpace* fes = positiveNodes->FESpace();
    const int numMeshElements = numElements();
    for(int elem=0; elem < numMeshElements; ++elem)
    {
      SpatialBoundingBox& bbox = m_elementBoundingBox[elem];

      // Add each dof of the element to the bbox
      // Note: positivity of Bernstein bases ensures that convex
      //       hull of element nodes contain entire element
      fes->GetElementDofs(elem, dofIndices);
      for(int i = 0; i< dofIndices.Size(); ++i)
      {
        int nIdx = dofIndices[i];

        SpacePoint pt;
        for(int j=0; j< NDIMS; ++j)
          pt[j] = (*positiveNodes)(fes->DofToVDof(nIdx,j));

        bbox.addPoint( pt );
      }

      // Slightly scale the bbox to account for numerical noise
      bbox.scale(bboxScaleFactor);

      m_meshBoundingBox.addBox( bbox );
    }

    /// Clean up -- deallocate grid function if necessary
    if(mustDeleteNodes)
    {
      delete positiveNodes;
      positiveNodes = AXOM_NULLPTR;
    }
  }

  /*!
   * Helper function to initialize the bounding boxes for a low order mfem mesh
   *
   * \param bboxScaleFactor Scale factor for expanding bounding boxes
   * \note This function can only be called before the class is initialized
   * \pre This function can only be called when m_isHighOrder is false
   * \pre bboxScaleFactor must be greater than or equal to 1.
   */
  void initializeBoundingBoxesLowOrder(double bboxScaleFactor)
  {
    SLIC_ASSERT( !m_isInitialized );
    SLIC_ASSERT( !m_isHighOrder );
    SLIC_ASSERT( bboxScaleFactor >= 1. );

    /// For each element, compute bounding box, and overall mesh bbox
    const int numMeshElements = numElements();
    for(int elem=0; elem < numMeshElements; ++elem)
    {
      SpatialBoundingBox& bbox = m_elementBoundingBox[elem];

      mfem::Element* elt = m_mesh->GetElement(elem);
      int* eltVerts = elt->GetVertices();
      for(int i = 0; i< elt->GetNVertices(); ++i)
      {
        int vIdx = eltVerts[i];
        bbox.addPoint( SpacePoint( m_mesh->GetVertex( vIdx ) ) );
      }

      // scale the bounding box to account for numerical noise
      bbox.scale(bboxScaleFactor);

      m_meshBoundingBox.addBox( bbox );
    }
  }

private:
  mfem::Mesh* m_mesh;
  bool m_isHighOrder;
  bool m_isInitialized;

  SpatialBoundingBox m_meshBoundingBox;
  std::vector<SpatialBoundingBox> m_elementBoundingBox;
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
 */
template < int NDIMS >
class PointInCell
{
public:

  typedef quest::ImplicitGrid<NDIMS> GridType;
  typedef typename GridType::SpacePoint SpacePoint;
  typedef typename GridType::SpatialBoundingBox SpatialBoundingBox;

  typedef int IndexType;

  /*!  Special value to indicate an unsuccessful query */
  static const IndexType NO_CELL;

  /*!
   * Construct a point in cell query structure over a given mfem mesh
   * \param mesh A pointer to an mfem mesh
   * \param resolution If the resolution is not provided, we use a heuristic to set the resolution
   */
  PointInCell(mfem::Mesh* mesh,
      const primal::Point<int, NDIMS>& resolution = primal::Point<int, NDIMS>::zero()) :
    m_meshWrapper(mesh), m_queryPointsMeshDumper(AXOM_NULLPTR), m_queryPathsMeshDumper(AXOM_NULLPTR)
  {
    int num_elem = m_meshWrapper.numElements();

    SpatialBoundingBox meshBB = m_meshWrapper.meshBoundingBox();
    SLIC_INFO("Mesh bounding box is: " << meshBB );

    // Slightly scale the bounding box
    double const EPS = 1e-8; 
    meshBB.scale(1. + EPS);
    m_grid.initialize(meshBB, resolution, num_elem);

    for(int i=0; i< num_elem; ++i)
    {
      SpatialBoundingBox elemBBox = m_meshWrapper.elementBoundingBox(i);
      m_grid.insert( elemBBox, i);
    }
  }

  ~PointInCell()
  {
    if(m_queryPathsMeshDumper != AXOM_NULLPTR)
    {
      delete m_queryPathsMeshDumper;
      m_queryPathsMeshDumper = AXOM_NULLPTR;
    }
  }

  void enableDebugMeshGeneration()
  {
    if(m_queryPathsMeshDumper == AXOM_NULLPTR)
    {
      m_queryPathsMeshDumper = new detail::PICQueryMeshDumper<NDIMS>;
    }
  }

  /*!
   * \brief Attempts to find the index of the mesh cell containing query point pt
   * \param[in] pt The query point
   * \param[out] isoparametric When pt is within a cell, this point returns
   *             the isoparametric coordinates of pt within the cell.
   * \return The index of the mesh cell containing pt, when one exists,
   *         otherwise returns the special value PointInCell::NO_CELL
   */
  IndexType locatePoint(const SpacePoint& pt, SpacePoint* isoparametric = AXOM_NULLPTR) const
  {
    typedef typename GridType::BitsetType BitsetType;

    IndexType containingCell = NO_CELL;
    static SpacePoint zero = SpacePoint::zero();
    SpacePoint& isopar = (isoparametric == AXOM_NULLPTR) ? zero : *isoparametric;

    BitsetType candidates = m_grid.getCandidates(pt);

    // SLIC_INFO("Candidates for space point " <<  pt);
    bool foundContainingCell = false;
    for(std::size_t eltIdx = candidates.find_first();
        !foundContainingCell && eltIdx != BitsetType::npos;
        eltIdx = candidates.find_next( eltIdx) )
    {
      // First check that pt is in bounding box of element
      if( m_meshWrapper.elementBoundingBox(eltIdx).contains(pt))
      {
        // if isopar is in the proper range
        if( getIsoparametricCoords(eltIdx, pt, &isopar) )
        {
          // then we have found the cellID
          foundContainingCell = true;
          containingCell = eltIdx;
        }
      }
    }

    return containingCell;
  }

  /*!
   * Returns the point in space corresponding to isoparameteric coordinates isopar
   * of the mesh element with index eltIdx
   */
  SpacePoint findInSpace(IndexType eltIdx, const SpacePoint& isopar)
  {
    SpacePoint pt;

    mfem::IsoparametricTransformation tr;
    m_meshWrapper.getMesh()->GetElementTransformation(eltIdx, &tr);

    mfem::IntegrationPoint ip;
    ip.Set(isopar.data(), NDIMS);

    mfem::Vector v(pt.data(), NDIMS);
    tr.Transform(ip, v);

    return pt;
  }

  /*!
   * Attempts to find isoparametric coordinates \a isopar of a point \a pt in space
   * with respect to a given element of the mesh (with index \a eltIdx)
   * \return True if \a pt is contained in the element, in which case
   *         the coordinates of \a isopar will be in the unit cube (of dimension NDIMS)
   */
  bool getIsoparametricCoords(IndexType eltIdx, const SpacePoint& pt, SpacePoint* isopar) const
  {
    const int refineOrder = 1;      // TODO: Make this a class variable

    mfem::IsoparametricTransformation tr;
    m_meshWrapper.getMesh()->GetElementTransformation(eltIdx, &tr);
    mfem::Vector ptSpace(const_cast<double*>(pt.data()), NDIMS);

    // TransformBack returns zero if the element is properly mapped
    mfem::IntegrationPoint ipRef;

    // Status codes: {0 -> successful; 1 -> pt was outside; 2-> did not converge}
    int err = tr.TransformBack(ptSpace, ipRef, refineOrder);
    ipRef.Get(isopar->data(), NDIMS);

    // Add query paths to debug mesh, if applicable
    if( m_queryPathsMeshDumper != AXOM_NULLPTR)
    {
      m_queryPathsMeshDumper->addQueryPoint(pt, *isopar, err, tr);
    }

    return (err == 0);
  }


  void printDebugMesh(const std::string& filename)
  {
    if( m_queryPathsMeshDumper != AXOM_NULLPTR)
    {
      m_queryPathsMeshDumper->outputDebugMesh(filename);
    }
  }

private:
  detail::PICMeshWrapper<NDIMS, mfem::Mesh> m_meshWrapper;
  GridType    m_grid;

  detail::PICQueryMeshDumper<NDIMS>* m_queryPointsMeshDumper;
  detail::PICQueryMeshDumper<NDIMS>* m_queryPathsMeshDumper;
};


template<int NDIMS>
const typename PointInCell<NDIMS>::IndexType PointInCell<NDIMS>::NO_CELL = -1;

} // end namespace quest
} // end namespace axom

#endif // QUEST_POINT_IN_CELL_HPP_
