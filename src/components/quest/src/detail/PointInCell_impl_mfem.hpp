/*
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#ifndef QUEST_POINT_IN_CELL_MFEM_IMPL_HPP_
#define QUEST_POINT_IN_CELL_MFEM_IMPL_HPP_

#include "axom/config.hpp"
#include "axom/Macros.hpp"

#include "slic/slic.hpp"

#include "primal/Point.hpp"
#include "primal/BoundingBox.hpp"

#ifdef AXOM_USE_MFEM
#  include "mfem.hpp"
#else
#  error PointInCell_mfem_impl depends on mfem.
#endif


namespace axom {
namespace quest {

/** Tag for mfem-based specialization of PointInCell query code */
struct quest_point_in_cell_mfem_tag {};


namespace detail {

template<typename quest_point_in_cell_mesh_tag>
struct PointInCellTraits;


template <typename quest_point_in_cell_mesh_tag>
class PointInCellMeshWrapper;


template<>
struct PointInCellTraits<quest_point_in_cell_mfem_tag>
{
public:
  typedef mfem::Mesh MeshType;
  typedef int IndexType;

  /*!  Special value to indicate an unsuccessful query */
  static const IndexType NO_CELL;
};

/** Initialize the NO_CELL static variable of PointInCellTraits for mfem */
const PointInCellTraits<quest_point_in_cell_mfem_tag>::IndexType
PointInCellTraits<quest_point_in_cell_mfem_tag>::NO_CELL = -1;



/*!
 * \brief Wraps MFEM mesh access functionality for the PointInCell class
 */
template <>
class PointInCellMeshWrapper< quest_point_in_cell_mfem_tag >
{
public:
  typedef int IndexType;
  typedef PointInCellMeshWrapper< quest_point_in_cell_mfem_tag > MeshWrapper;

  PointInCellMeshWrapper(mfem::Mesh* mesh) : m_mesh(mesh)
  {
    // Some sanity checks
    SLIC_ASSERT( m_mesh != AXOM_NULLPTR);

    m_isHighOrder = (m_mesh->GetNodalFESpace() != AXOM_NULLPTR) && (m_mesh->GetNE() > 0);

/*
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
*/
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
  static mfem::FiniteElementCollection* getCorrespondingPositiveFEC(
      const mfem::FiniteElementCollection* fec,
      int order,
      int dim,
      int mapType)
  {
    //   NOTE(KW) Should this use a map<string, string> instead ?
    //         It could then use FECollection::New(string)

   // If fec is already positive, return it?
    SLIC_CHECK_MSG( !isPositiveBasis(fec),
        "This function is only meant to be called on non-positive finite element collection"
        );

    // Attempt to find the corresponding positive H1 fec
    if(dynamic_cast<const mfem::H1_FECollection*>(fec))
    {
      return new mfem::H1_FECollection(order, dim, mfem::BasisType::Positive);
    }

    // Attempt to find the corresponding positive L2 fec
    if(dynamic_cast<const mfem::L2_FECollection*>(fec))
    {
      // should we throw a not supported error here?
      return new mfem::L2_FECollection(order, dim, mfem::BasisType::Positive, mapType);
    }

    // Attempt to find the corresponding quadratic or cubic fec
    if(dynamic_cast<const mfem::QuadraticFECollection*>(fec) ||
       dynamic_cast<const mfem::CubicFECollection*>(fec) )
    {
      SLIC_ASSERT( order == 2 || order == 3);
      return new mfem::H1_FECollection(order, dim, mfem::BasisType::Positive);
    }

    // Give up -- return NULL
    return AXOM_NULLPTR;
  }


  /*! Get the number of elements in the mesh */
  int numElements() const { return m_mesh->GetNE(); }

  int meshDimension() const { return m_mesh->Dimension(); }

  /*! Get a pointer to the mesh */
  mfem::Mesh* getMesh() const { return m_mesh; }

  template<int NDIMS>
  void computeBoundingBoxes(double bboxScaleFactor,
      std::vector< axom::primal::BoundingBox<double, NDIMS> >& eltBBoxes,
      axom::primal::BoundingBox<double, NDIMS>& meshBBox)
  {
    SLIC_ASSERT( static_cast<int>(eltBBoxes.size()) >= numElements() );

    if(m_isHighOrder)
    {
      computeHighOrderBoundingBoxes(bboxScaleFactor, eltBBoxes, meshBBox);
    }
    else
    {
      computeLowOrderBoundingBoxes(bboxScaleFactor, eltBBoxes, meshBBox);
    }
  }

  void findInSpace(IndexType eltIdx, const double* isopar, double* pt)
  {
    const int dim = meshDimension();

    mfem::IsoparametricTransformation tr;
    m_mesh->GetElementTransformation(eltIdx, &tr);

    mfem::IntegrationPoint ip;
    ip.Set(isopar, dim);

    mfem::Vector v(pt, dim);
    tr.Transform(ip, v);
  }

  /*!
   * Attempts to find isoparametric coordinates \a isopar of a point \a pt in space
   * with respect to a given element of the mesh (with index \a eltIdx)
   * \return True if \a pt is contained in the element, in which case
   *         the coordinates of \a isopar will be in the unit cube (of dimension NDIMS)
   */
  bool getIsoparametricCoords(IndexType eltIdx, const double* pt, double* isopar)
  {
    const int dim = meshDimension();
    const int refineOrder = 1;      // TODO: Make this a class variable

    mfem::IsoparametricTransformation tr;
    m_mesh->GetElementTransformation(eltIdx, &tr);
    mfem::Vector ptSpace(const_cast<double*>(pt), dim);

    // TransformBack returns zero if the element is properly mapped
    mfem::IntegrationPoint ipRef;

    // Status codes: {0 -> successful; 1 -> pt was outside; 2-> did not converge}
    int err = tr.TransformBack(ptSpace, ipRef, refineOrder);
    ipRef.Get(isopar, dim);

//    // Add query paths to debug mesh, if applicable
//    if( m_queryPathsMeshDumper != AXOM_NULLPTR)
//    {
//      m_queryPathsMeshDumper->addQueryPoint(pt, *isopar, err, tr);
//    }

    return (err == 0);
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
  template<int NDIMS>
  void computeHighOrderBoundingBoxes(double bboxScaleFactor,
      std::vector< axom::primal::BoundingBox<double, NDIMS> >& eltBBoxes,
      axom::primal::BoundingBox<double, NDIMS>& meshBBox)
  {
    typedef axom::primal::Point<double, NDIMS> SpacePoint;
    typedef axom::primal::BoundingBox<double, NDIMS> SpatialBoundingBox;

    // Sanity checks
    SLIC_ASSERT( m_isHighOrder );
    SLIC_ASSERT(bboxScaleFactor >= 1.);

    /// Generate (or access existing) positive (Bernstein) nodal grid function
    const mfem::FiniteElementSpace* nodalFESpace = m_mesh->GetNodalFESpace();

    mfem::GridFunction* positiveNodes = AXOM_NULLPTR;
    bool mustDeleteNodes = false;

    const mfem::FiniteElementCollection* nodalFEColl = nodalFESpace->FEColl();

    // Check if grid function is positive, if not create positive grid function
    if( MeshWrapper::isPositiveBasis( nodalFEColl ) )
    {
      positiveNodes = m_mesh->GetNodes();
    }
    else
    {
      // Assume that all elements of the mesh have the same order and geom type
      int order = nodalFESpace->GetOrder(0);
      int dim = meshDimension();
      int GeomType = m_mesh->GetElementBaseGeometry(0);
      int mapType = (nodalFEColl != AXOM_NULLPTR)
          ? nodalFEColl->FiniteElementForGeometry(GeomType)->GetMapType()
          : static_cast<int>(mfem::FiniteElement::VALUE);

      mfem::FiniteElementCollection* posFEColl =
          MeshWrapper::getCorrespondingPositiveFEC(nodalFEColl, order, dim, mapType);

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
      SpatialBoundingBox& bbox = eltBBoxes[elem];

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

      meshBBox.addBox( bbox );
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
  template<int NDIMS>
  void computeLowOrderBoundingBoxes(double bboxScaleFactor,
      std::vector< axom::primal::BoundingBox<double, NDIMS> >& eltBBoxes,
      axom::primal::BoundingBox<double, NDIMS>& meshBBox)
  {
    typedef axom::primal::Point<double, NDIMS> SpacePoint;
    typedef axom::primal::BoundingBox<double, NDIMS> SpatialBoundingBox;

    SLIC_ASSERT( !m_isHighOrder );
    SLIC_ASSERT( bboxScaleFactor >= 1. );

    /// For each element, compute bounding box, and overall mesh bbox
    const int numMeshElements = numElements();
    for(int elem=0; elem < numMeshElements; ++elem)
    {
      SpatialBoundingBox& bbox = eltBBoxes[elem];

      mfem::Element* elt = m_mesh->GetElement(elem);
      int* eltVerts = elt->GetVertices();
      for(int i = 0; i< elt->GetNVertices(); ++i)
      {
        int vIdx = eltVerts[i];
        bbox.addPoint( SpacePoint( m_mesh->GetVertex( vIdx ) ) );
      }

      // scale the bounding box to account for numerical noise
      bbox.scale(bboxScaleFactor);

      meshBBox.addBox( bbox );
    }
  }

private:
  mfem::Mesh* m_mesh;
  bool m_isHighOrder;
};




} // end namespace detail

} // end namespace quest
} // end namespace axom



#endif // QUEST_POINT_IN_CELL_MFEM_IMPL_HPP_
