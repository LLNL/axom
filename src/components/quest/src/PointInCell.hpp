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
#  error Mfem is a required dependency of quest::PointInCell.
#endif

#include "mfem.hpp"

namespace axom {
namespace quest {


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

  PICMeshWrapper(mfem::Mesh* mesh) :
    m_mesh(mesh), m_bernsteinNodes(AXOM_NULLPTR), m_mustDeleteNodes(false)
  {
    // Some sanity checks
    SLIC_ASSERT( m_mesh != AXOM_NULLPTR);
    SLIC_ASSERT( m_mesh->Dimension() == NDIMS);

    // Generate (or access existing) positive (Bernstein) nodal grid function
    const mfem::FiniteElementSpace* nodalFESpace = m_mesh->GetNodalFESpace();
    m_isHighOrder = (nodalFESpace != AXOM_NULLPTR) && (m_mesh->GetNE() > 0);

    if( m_isHighOrder )
    {
      // Assume that all elements of the mesh have the same order
      const mfem::FiniteElementCollection* nodalFEColl = nodalFESpace->FEColl();
      int order = nodalFESpace->GetOrder(0);

      bool isPositive = PICMeshWrapper::isPositiveBasis( nodalFEColl );
      if(isPositive)
      {
        m_bernsteinNodes = m_mesh->GetNodes();
      }
      else
      {
        int GeomType = m_mesh->GetElementBaseGeometry(0);
        int mapType = (nodalFEColl != AXOM_NULLPTR)
            ? nodalFEColl->FiniteElementForGeometry(GeomType)->GetMapType()
            : static_cast<int>(mfem::FiniteElement::VALUE);

        mfem::FiniteElementCollection* posFEColl =
            PICMeshWrapper::getCorrespondingPositiveFEC(nodalFEColl, order, mapType);

        if(posFEColl != AXOM_NULLPTR)
        {
          // Create a positive (Bernstein) grid function for the nodes
          mfem::FiniteElementSpace* posFESpace =
              new mfem::FiniteElementSpace(m_mesh, posFEColl, NDIMS);
          m_bernsteinNodes = new mfem::GridFunction(posFESpace);
          m_bernsteinNodes->MakeOwner(posFEColl);
          m_mustDeleteNodes = true;

          // Project the nodal grid function onto this
          m_bernsteinNodes->ProjectGridFunction(*(m_mesh->GetNodes()));
        }
      }

      SLIC_INFO("Mesh nodes fec -- " << nodalFEColl->Name()
          << " with ordering " << nodalFESpace->GetOrdering()
          << "\n\t -- Positive nodes are fec -- "
          << m_bernsteinNodes->FESpace()->FEColl()->Name()
          << " with ordering " << m_bernsteinNodes->FESpace()->GetOrdering() );


    }
  }

  ~PICMeshWrapper()
  {
    if(m_mustDeleteNodes)
    {
      delete m_bernsteinNodes;
      m_bernsteinNodes = AXOM_NULLPTR;
    }
  }


  /*! Predicate to check if the given basis is positive (Bernstein) */
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
   * \brief Utility function to get a positive (Bernstein) basis corresponding to the given
   *        FiniteElementCollection.
   */
  static mfem::FiniteElementCollection* getCorrespondingPositiveFEC(const mfem::FiniteElementCollection* fec, int order, int mapType)
  {
    //   NOTE(KW) Should this use a map<string, string> instead ?
    //         It could then use FECollection::New(string)

   // If fec is already positive, return it
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

  /*! Returns the vertex (low order) or node (high order) at the given index */
  SpacePoint getNode(int idx) const
  {
    if(m_isHighOrder)
    {
      SpacePoint pt;

      mfem::FiniteElementSpace* fes = m_bernsteinNodes->FESpace();
      for(int j=0; j< NDIMS; ++j)
        pt[j] = (*m_bernsteinNodes)(fes->DofToVDof(idx,j));

      return pt;
    }
    else
    {
      return SpacePoint(m_mesh->GetVertex(idx) );
    }

  }

  /*! Compute the bounding box of the mesh */
  SpatialBoundingBox meshBoundingBox() const
  {
    SpatialBoundingBox bbox;

    int numNodes = m_isHighOrder
        ? m_bernsteinNodes->FESpace()->GetNDofs()
        : m_mesh->GetNV();

    // Add each vertex (low order) or node (high order) to the bbox
    for(int i=0; i< numNodes; ++i)
    {
      bbox.addPoint( getNode(i) );
    }

    return bbox;
  }

  /*! Compute the bounding box of mesh element with index idx */
  SpatialBoundingBox elementBoundingBox(IndexType idx) const
  {
    SpatialBoundingBox bbox;

    if(m_isHighOrder)
    {
      mfem::Array<int> dofIndices;

      m_bernsteinNodes->FESpace()->GetElementDofs(idx, dofIndices);
      for(int i = 0; i< dofIndices.Size(); ++i)
      {
        int nIdx = dofIndices[i];
        bbox.addPoint( getNode(nIdx) );
      }
    }
    else
    {
      mfem::Element* elt = m_mesh->GetElement(idx);
      int* eltVerts = elt->GetVertices();
      for(int i = 0; i< elt->GetNVertices(); ++i)
      {
        int vIdx = eltVerts[i];
        bbox.addPoint( getNode( vIdx ) );
      }
    }

    double const EPS = 1e-8;  // TODO -- convert to class member ?
    bbox.scale(1. + EPS);

    return bbox;
  }


private:
  mfem::Mesh* m_mesh;
  mfem::GridFunction* m_bernsteinNodes;
  bool m_isHighOrder;
  bool m_mustDeleteNodes;
};



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
  PointInCell(mfem::Mesh* mesh, const primal::Point<int, NDIMS>& resolution = primal::Point<int, NDIMS>::zero()) :
    m_meshWrapper(mesh), m_transformBackDebugMesh(2), m_numQueries(0)
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
   * Attempts to find the isoparametric coordinates isopar of a point in space pt,
   * with respect to a given element of the mesh (with index eltIdx)
   * \return True if the pt is contained in the element, in which case
   *         the coordinates of isopar will be in the unit cube (of dimension NDIMS)
   */
  bool getIsoparametricCoords(IndexType eltIdx, const SpacePoint& pt, SpacePoint* isopar) const
  {
    mfem::IsoparametricTransformation tr;
    m_meshWrapper.getMesh()->GetElementTransformation(eltIdx, &tr);
    mfem::Vector ptSpace(const_cast<double*>(pt.data()), NDIMS);

    // TransformBack returns zero if the element is properly mapped
    mfem::IntegrationPoint ipRef;

    int err = tr.TransformBack(ptSpace, ipRef);

    int meshStartIndex = m_transformBackDebugMesh.getMeshNumberOfNodes();

    {
      //SLIC_INFO("Used " << tr.tbDebug.numIters << " iterations." );

      // Add vertices, including the last iteration's results
      mfem::Vector y;
      tr.Transform(ipRef, y);
      tr.tbDebug.points.SetCol( tr.tbDebug.numIters, y );

      for(int i=0; i< tr.tbDebug.numIters; ++i)
      {
        double x = tr.tbDebug.points(0, i);
        double y = tr.tbDebug.points(1, i);
        m_transformBackDebugMesh.insertNode(x,y);
      }

      for(int i=0; i< tr.tbDebug.numIters-1; ++i)
      {
        int currVertIdx = meshStartIndex + i;
        int indices[2] = { currVertIdx, currVertIdx+1 };
        m_transformBackDebugMesh.insertCell(indices, MINT_SEGMENT, 2);
        m_tb_query_index.push_back( m_numQueries );
        m_tb_query_return_code.push_back( err );
        m_tb_query_num_iters.push_back( tr.tbDebug.numIters );
        m_tb_query_iter.push_back(i);
      }

      m_numQueries++;
    }

    switch(err)
    {
    case 0: // success
      ipRef.Get(isopar->data(), NDIMS);

      return true;
    case 1: // out-of-bounds, early return
//      SLIC_INFO("Early return: Point " << pt << " was not found in element " << eltIdx);
      return false;
    case 2: // max-iter
//      SLIC_INFO("Point " << pt << " was not found in element " << eltIdx << " after 16 iterations.");
      return false;
    }

    return false;
  }

  void printDebugMesh(const std::string& filename)
  {
    mint::FieldData* CD = m_transformBackDebugMesh.getCellFieldData();

    int size = m_tb_query_index.size();

    if(size == 0)
    {
        SLIC_INFO("No mesh data.");
        return;
    }

    {
      std::string name = "query_index";
      CD->addField( new mint::FieldVariable< int >(name, size) );
      int* fld = CD->getField( name )->getIntPtr();

      for(int i=0; i < size; ++i)
      {
        fld[i] = m_tb_query_index[i];
      }
    }
    {
      std::string name = "query_status";
      CD->addField( new mint::FieldVariable< int >(name, size) );
      int* fld = CD->getField( name )->getIntPtr();

      for(int i=0; i < size; ++i)
      {
        fld[i] = m_tb_query_return_code[i];
      }
    }
    {
      std::string name = "query_num_iters";
      CD->addField( new mint::FieldVariable< int >(name, size) );
      int* fld = CD->getField( name )->getIntPtr();

      for(int i=0; i < size; ++i)
      {
        fld[i] = m_tb_query_num_iters[i];
      }
    }
    {
      std::string name = "query_iter";
      CD->addField( new mint::FieldVariable< int >(name, size) );
      int* fld = CD->getField( name )->getIntPtr();

      for(int i=0; i < size; ++i)
      {
        fld[i] = m_tb_query_iter[i];
      }
    }

    mint::write_vtk(&m_transformBackDebugMesh, filename);
  }

private:
  PICMeshWrapper<NDIMS, mfem::Mesh> m_meshWrapper;
  GridType    m_grid;

  mutable mint::UnstructuredMesh<MINT_SEGMENT> m_transformBackDebugMesh;
  mutable std::vector<int> m_tb_query_index;
  mutable std::vector<int> m_tb_query_return_code;
  mutable std::vector<int> m_tb_query_num_iters;
  mutable std::vector<int> m_tb_query_iter;
  mutable int m_numQueries;



};

template<int NDIMS>
const typename PointInCell<NDIMS>::IndexType PointInCell<NDIMS>::NO_CELL = -1;

} // end namespace quest
} // end namespace axom

#endif // QUEST_POINT_IN_CELL_HPP_
