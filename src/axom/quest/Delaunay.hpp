
#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"
#include "axom/spin.hpp"

#include "axom/fmt.hpp"

#include <list>
#include <vector>
#include <set>
#include <cstdlib>

namespace axom
{
namespace quest
{
/**
 * \brief A class for incremental generation of a 2D or 3D Delaunay triangulation
 *
 * Construct a Delaunay triangulation incrementally by inserting points one by one.
 * A bounding box of the points needs to be defined first via \a initializeBoundary(...)
 */
template <int DIM = 2>
class Delaunay
{
public:
  AXOM_STATIC_ASSERT_MSG(DIM == 2 || DIM == 3,
                         "The template parameter DIM can only be 2 or 3. ");

  using DataType = double;

  using PointType = primal::Point<DataType, DIM>;
  using ElementType =
    typename std::conditional<DIM == 2,
                              primal::Triangle<DataType, 2>,
                              primal::Tetrahedron<DataType, 3>>::type;
  using BaryCoordType = primal::Point<DataType, DIM + 1>;
  using BoundingBox = primal::BoundingBox<DataType, DIM>;

  using IAMeshType = slam::IAMesh<DIM, DIM, PointType>;
  using IndexType = typename IAMeshType::IndexType;
  using IndexArray = typename IAMeshType::IndexArray;

  using IndexPairType = std::pair<IndexType, IndexType>;

  static constexpr int VERT_PER_ELEMENT = DIM + 1;
  static constexpr IndexType INVALID_INDEX = -1;

private:
  using ModularFaceIndex =
    slam::ModularInt<slam::policies::CompileTimeSize<IndexType, VERT_PER_ELEMENT>>;

private:
  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;

public:
  /**
   * \brief Default constructor
   * \note User need to call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay()
    : m_has_boundary(false)
    , m_num_removed_elements_since_last_compact(0)
  { }

  /**
   * \brief Defines the boundary of the triangulation.
   * \details subsequent points added to the triangulation must not be outside of this boundary.
   */
  void initializeBoundary(const BoundingBox& bb)
  {
    std::vector<DataType> points;
    IndexArray elem;

    generateInitialMesh(points, elem, bb);

    m_mesh = IAMeshType(points, elem);

    m_bounding_box = bb;
    m_has_boundary = true;
  }

  /**
   * \brief Adds a new point and locally re-triangulates the mesh to ensure that it stays Delaunay
   *
   * This function will traverse the mesh to find the element that contains
   * this point, creates the Delaunay cavity, which takes out all the elements
   * that contains the point in its sphere, and fill it with a Delaunay ball.
   *
   * \pre The current mesh must already be Delaunay.
   */
  void insertPoint(const PointType& new_pt)
  {
    //Make sure initializeBoundary(...) is called first
    SLIC_ASSERT_MSG(
      m_has_boundary,
      "Error: Need a predefined boundary box prior to adding points.");

    //Make sure the new point is inside the boundary box
    SLIC_ASSERT_MSG(m_bounding_box.contains(new_pt),
                    "Error: new point is outside of the boundary box.");

    // Find the mesh element containing the insertion point
    IndexType element_i = findContainingElement(new_pt);

    if(element_i == INVALID_INDEX)
    {
      SLIC_WARNING(
        fmt::format("Could not insert point {} into Delaunay triangulation: "
                    "Element containing that point was not found",
                    new_pt));
      return;
    }

    // Run the insertion operation by finding invalidated elements around the point (the "cavity")
    // and replacing them with new valid elements (the Delaunay "ball")
    InsertionHelper insertionHelper(m_mesh);
    insertionHelper.findCavityElements(new_pt, element_i);
    insertionHelper.createCavity();
    IndexType new_pt_i = m_mesh.addVertex(new_pt);
    insertionHelper.delaunayBall(new_pt_i);

    m_num_removed_elements_since_last_compact +=
      insertionHelper.numRemovedElements();

    // Compact the mesh if there are too many removed elements
    if(shouldCompactMesh())
    {
      this->compactMesh();
    }
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 2, ElementType>::type getElement(
    int element_index) const
  {
    const auto verts = m_mesh.ev_rel[element_index];
    const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
    const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
    const PointType& p2 = m_mesh.getVertexPoint(verts[2]);

    return ElementType(p0, p1, p2);
  }

  template <int TDIM = DIM>
  typename std::enable_if<TDIM == 3, ElementType>::type getElement(
    int element_index) const
  {
    const auto verts = m_mesh.ev_rel[element_index];
    const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
    const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
    const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
    const PointType& p3 = m_mesh.getVertexPoint(verts[3]);

    return ElementType(p0, p1, p2, p3);
  }

  /**
   * \brief Prints out mesh details, for debugging purpose.
   */
  void printMesh() { m_mesh.print_all(); }

  /**
   * \brief Write the m_mesh to a legacy VTK file
   *
   * \param filename The name of the file to write to,
   * \note The suffix ".vtk" will be appended to the provided filename
   * \details This function uses mint to write the m_mesh to VTK format.
   */
  void writeToVTKFile(const std::string& filename)
  {
    const auto CELL_TYPE = DIM == 2 ? mint::TRIANGLE : mint::TET;
    mint::UnstructuredMesh<mint::SINGLE_SHAPE> mint_mesh(DIM, CELL_TYPE);

    this->compactMesh();

    for(auto v : m_mesh.vertex_set.positions())
    {
      mint_mesh.appendNodes(m_mesh.getVertexPoint(v).data(), 1);
    }

    for(auto e : m_mesh.element_set.positions())
    {
      mint_mesh.appendCell(&(m_mesh.ev_rel[e][0]), CELL_TYPE);
    }

    mint::write_vtk(&mint_mesh, filename);
  }

  /**
   * \brief Removes the vertices that defines the boundary of the mesh,
   * and the elements attached to them.
   *
   * \details After this function is called, no more points can be added to the m_mesh.
   */
  void removeBoundary()
  {
    if(m_has_boundary)
    {
      //remove the boundary box, which will be the first 4 points for triangles, first 8 for tetrahedron
      const int num_boundary_pts = 1 << DIM;

      //Collect a list of elements to remove first, because
      //the list may be incomplete if generated during the removal.
      IndexArray elements_to_remove;
      for(int v = 0; v < num_boundary_pts; ++v)
      {
        IndexArray elems = m_mesh.getElementsWithVertex(v);
        elements_to_remove.insert(elements_to_remove.end(),
                                  elems.begin(),
                                  elems.end());
      }
      for(auto e : elements_to_remove)
      {
        if(m_mesh.isValidElementEntry(e))
        {
          m_mesh.removeElement(e);
        }
      }

      for(int v = 0; v < num_boundary_pts; ++v)
      {
        m_mesh.removeVertex(v);
      }

      this->compactMesh();
      m_has_boundary = false;
    }
  }

  /// \brief Get the IA mesh data pointer
  const IAMeshType* getMeshData() const { return &m_mesh; }

  /**
   * \brief Checks that the underlying mesh is a valid Delaunay triangulation of the point set
   *
   * A Delaunay triangulation is valid when none of the vertices are inside the circumspheres 
   * of any of the elements of the mesh
   */
  bool isValid(bool verboseOutput = false) const
  {
    // We use an ImplicitGrid spatial index to find the candidate elements
    // whose circumsphere might contain the vertices of the mesh

    using ImplicitGridType = spin::ImplicitGrid<DIM>;
    using NumericArrayType = primal::NumericArray<DataType, DIM>;
    using BitsetType = typename ImplicitGridType::BitsetType;

    bool valid = true;

    std::vector<std::pair<IndexType, IndexType>> invalidEntries;

    ImplicitGridType grid(m_bounding_box, nullptr, m_mesh.element_set.size());

    // Add (bounding boxes of) element circumspheres to implicit grid
    for(auto element_idx : m_mesh.element_set)
    {
      if(m_mesh.isValidElementEntry(element_idx))
      {
        const auto sphere = this->getElement(element_idx).circumsphere();
        const auto center = NumericArrayType(sphere.getCenter());
        const auto offset = NumericArrayType(sphere.getRadius());

        BoundingBox bb;
        bb.addPoint(PointType(center - offset));
        bb.addPoint(PointType(center + offset));

        grid.insert(bb, element_idx);  // insert valid entries into grid
      }
    }

    // for each vertex -- check in_sphere condition for candidate element
    for(auto vertex_idx : m_mesh.vertex_set)
    {
      const auto& vertex = m_mesh.getVertexPoint(vertex_idx);

      auto candidates = grid.getCandidates(vertex);
      for(auto element_idx = candidates.find_first();  //
          element_idx != BitsetType::npos;
          element_idx = candidates.find_next(element_idx))
      {
        // skip if element at this index is not valid
        if(!m_mesh.isValidElementEntry(element_idx))
        {
          continue;
        }

        // also skip if this is a vertex of the element
        if(slam::is_subset<IndexType>(vertex_idx,
                                      m_mesh.getVerticesInElement(element_idx)))
        {
          continue;
        }

        // check insphere condition
        if(primal::in_sphere(vertex, getElement(element_idx)))
        {
          valid = false;

          if(verboseOutput)
          {
            invalidEntries.push_back(std::make_pair(vertex_idx, element_idx));
          }
        }
      }
    }

    if(verboseOutput)
    {
      if(valid)
      {
        SLIC_INFO("Delaunay complex was valid");
      }
      else
      {
        fmt::memory_buffer out;
        for(const auto& pr : invalidEntries)
        {
          const auto vertex_idx = pr.first;
          const auto element_idx = pr.second;
          const auto& pos = m_mesh.getVertexPoint(vertex_idx);
          const auto element = this->getElement(element_idx);
          const auto circumsphere = element.circumsphere();
          fmt::format_to(out,
                         "\n\tVertex {} @ {}"
                         "\n\tElement {}: {} w/ circumsphere: {}"
                         "\n\tDistance to circumcenter: {}",
                         vertex_idx,
                         pos,
                         element_idx,
                         element,
                         circumsphere,
                         std::sqrt(primal::squared_distance(
                           pos,
                           PointType(circumsphere.getCenter()))));
        }

        SLIC_INFO(
          fmt::format("Delaunay complex was NOT valid. There were {} "
                      "vertices in the circumsphere of an element. {}",
                      invalidEntries.size(),
                      fmt::to_string(out)));
      }
    }

    return valid;
  }

private:
  /// \brief Find the index of the element that contains the query point, or the element closest to the point.
  IndexType findContainingElement(const PointType& query_pt)
  {
    if(m_mesh.isEmpty())
    {
      SLIC_ERROR(
        "Attempting to insert point into empty Delaunay triangulation."
        "Delaunay::initializeBoundary() needs to be called first");
      return INVALID_INDEX;
    }
    if(!m_bounding_box.contains(query_pt))
    {
      SLIC_WARNING(
        "Attempting to locate element at location outside valid domain");
      return INVALID_INDEX;
    }

    //find the last valid element to use as starting element
    IndexType element_i = m_mesh.getValidElementIndex();

    while(1)
    {
      const BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);

      //Find the index of the most negative barycentric coord
      //Use modular index since it could wrap around to 0
      ModularFaceIndex modular_idx(bary_coord.array().argMin());

      if(bary_coord[modular_idx] >= 0)  // inside if smallest bary coord positive
      {
        return element_i;
      }

      // else, move to that neighbor
      element_i = m_mesh.ee_rel[element_i][modular_idx + 1];

      // Either there is a hole in the m_mesh, or the point is outside of the m_mesh.
      // Logically, this should never happen.
      if(!m_mesh.isValidElementEntry(element_i))
      {
        SLIC_WARNING(fmt::format(
          "Entered invalid element in "
          "Delaunay::findContainingElement(). Underlying mesh {} valid",
          m_mesh.isValid() ? "is" : "is not"));
        return INVALID_INDEX;
      }
    }
  }

  /// \brief Predicate for when to compact internal mesh data structures after removing elements
  bool shouldCompactMesh() const
  {
    // Note: This auto-compacting feature is hard coded.
    // It may be good to let user have control of this option in the future.
    return m_num_removed_elements_since_last_compact > 512 &&
      (m_num_removed_elements_since_last_compact > .2 * m_mesh.element_set.size());
  }

  /// \brief Compacts the underlying mesh
  void compactMesh()
  {
    m_mesh.compact();
    m_num_removed_elements_since_last_compact = 0;
  }

  /**
   * \brief Helper function to fill the array with the initial mesh.
   * \details create a rectangle for 2D, cube for 3D, and fill the array with the mesh data.
   */
  void generateInitialMesh(std::vector<DataType>& points,
                           std::vector<IndexType>& elem,
                           const BoundingBox& bb);

  /**
   * \brief helper function to retrieve the barycentric coordinate of the query point in the element
   */
  BaryCoordType getBaryCoords(IndexType element_idx, const PointType& q_pt) const;

private:
  /**
   * \brief Helper class to locally insert a new point into a Delaunay complex while keeping the mesh Delaunay
   */
  struct InsertionHelper
  {
  public:
    InsertionHelper(IAMeshType& mesh)
      : m_mesh(mesh)
      , facet_set(0)
      , fv_rel(&facet_set, &m_mesh.vertex_set)
      , fc_rel(&facet_set, &m_mesh.element_set)
      , cavity_elems(0)
      , inserted_elems(0)
    { }

    /**
   * \brief Find the list of elements whose circumspheres contain the query point
   *
   * \details This function starts from an element \a element_i, and searches through
   * neighboring elements for a list of element indices whose circumspheres contains the query point.
   * It uses a recursive helper function \a findCavityElementsRec()
   *
   * \param query_pt the query point
   * \param element_i the element to start the search at
   */
    void findCavityElements(const PointType& query_pt, IndexType element_i)
    {
      m_checked_element_set.insert(element_i);
      findCavityElementsRec(query_pt, element_i);

      SLIC_ASSERT_MSG(!cavity_elems.empty(),
                      "Error: New point is not contained in the mesh");
      SLIC_ASSERT(!facet_set.empty());
    }

    /**
   * \brief recursive function to find cavity elements given a point to be added
   * \details Check if the point is in the circumsphere element, if so,
   * recursively call the neighboring elements.
   */
    bool findCavityElementsRec(const PointType& query_pt, IndexType element_idx)
    {
      const bool is_in_circle = isPointInSphere(query_pt, element_idx);

      // Base case: Element element_idx was valid (i.e. query_pt is outside its circumcircle)
      if(!is_in_circle)
      {
        return true;
      }

      // Recursive case: Add element to list and check its neighbors
      cavity_elems.insert(element_idx);

      //check for each faces; m_checked_element_set ensures we only check each element once
      for(auto nbr = m_mesh.ee_rel.begin(element_idx);
          nbr != m_mesh.ee_rel.end(element_idx);
          ++nbr)
      {
        //The latter case is a checked element that is not a cavity element
        if(!m_mesh.isValidElementEntry(*nbr) ||
           (m_checked_element_set.insert(*nbr).second
              ? findCavityElementsRec(query_pt, *nbr)
              : cavity_elems.findIndex(*nbr) == ElementSet::INVALID_ENTRY))
        {
          IndexArray vlist = m_mesh.getElementFace(element_idx, nbr.index());

          //For tetrahedron, if the element face is odd, reverse vertex order
          if(DIM == 3 && nbr.index() % 2 == 1)
          {
            axom::utilities::swap(vlist[1], vlist[2]);
          }

          auto fIdx = facet_set.insert();
          for(int i = 0; i < VERTS_PER_FACET; i++)
          {
            fv_rel.insert(fIdx, vlist[i]);
          }

          fc_rel.insert(fIdx, *nbr);
        }
      }
      return false;
    }

    /**
    * \brief Remove the elements in the Delaunay cavity
    */
    void createCavity()
    {
      for(auto elem : cavity_elems)
      {
        m_mesh.removeElement(elem);
      }
    }

    /// \brief Fill in the Delaunay cavity with new elements containing the insertion point
    void delaunayBall(IndexType new_pt_i)
    {
      const int numFaces = facet_set.size();

      for(int i = 0; i < numFaces; ++i)
      {
        // Create a new element from the face and the inserted point
        IndexType vlist[VERT_PER_ELEMENT];
        for(int d = 0; d < VERTS_PER_FACET; ++d)
        {
          vlist[d] = fv_rel[i][d];
        }
        vlist[VERTS_PER_FACET] = new_pt_i;

        // set all neighbors to nID; they'll be fixed in the fixVertexNeighborhood function below
        IndexType neighbors[VERT_PER_ELEMENT];
        const auto nID = fc_rel[i][0];
        for(int d = 0; d < VERTS_PER_FACET; ++d)
        {
          neighbors[d] = nID;
        }

        IndexType new_el = m_mesh.addElement(vlist, neighbors);
        inserted_elems.insert(new_el);
      }

      // Fix neighborhood around the new point
      m_mesh.fixVertexNeighborhood(new_pt_i, inserted_elems.data());
    }

    /// \brief Returns the number of elements removed during this insertion
    int numRemovedElements() const { return cavity_elems.size(); }

    /// \brief Helper function returns true if the query point is in the sphere formed by the element vertices
    bool isPointInSphere(const PointType& query_pt, IndexType element_idx) const;

  public:
    // we create a surface mesh
    // sets: vertex, facet
    using PositionType = typename IAMeshType::PositionType;
    using ElementType = typename IAMeshType::ElementType;

    using ElementSet = typename IAMeshType::ElementSet;
    using VertexSet = typename IAMeshType::VertexSet;
    using FacetSet = slam::DynamicSet<PositionType, ElementType>;

    // relations: facet->vertex, facet->cell
    static constexpr int VERTS_PER_FACET = IAMeshType::VERTS_PER_ELEM - 1;
    using FacetBoundaryRelation =
      typename IAMeshType::template IADynamicConstantRelation<VERTS_PER_FACET>;
    using FacetCoboundaryRelation =
      typename IAMeshType::template IADynamicConstantRelation<1>;

  public:
    IAMeshType& m_mesh;

    FacetSet facet_set;
    FacetBoundaryRelation fv_rel;
    FacetCoboundaryRelation fc_rel;

    ElementSet cavity_elems;
    ElementSet inserted_elems;

    std::set<IndexType> m_checked_element_set;
  };
};

//--------------------------------------------------------------------------------
// Below are 2D and 3D specializations for methods in the Delaunay class
//--------------------------------------------------------------------------------

// 2D specialization for generateInitialMesh(...)
template <>
void Delaunay<2>::generateInitialMesh(std::vector<DataType>& points,
                                      std::vector<IndexType>& elem,
                                      const BoundingBox& bb)
{
  //Set up the initial IA mesh of 2 triangles forming a rectangle

  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  // clang-format off
  std::vector<DataType> pt { mins[0], mins[1], 
                             mins[0], maxs[1], 
                             maxs[0], mins[1], 
                             maxs[0], maxs[1] };

  std::vector<IndexType> el { 0, 2, 1, 
                              3, 1, 2 };
  // clang-format on

  points.swap(pt);
  elem.swap(el);
}

// 3D specialization for generateInitialMesh(...)
template <>
void Delaunay<3>::generateInitialMesh(std::vector<DataType>& points,
                                      std::vector<IndexType>& elem,
                                      const BoundingBox& bb)
{
  //Set up the initial IA mesh of 6 tetrahedrons forming a cube
  const PointType& mins = bb.getMin();
  const PointType& maxs = bb.getMax();

  // clang-format off
  std::vector<DataType> pt { mins[0], mins[1], mins[2], 
                             mins[0], mins[1], maxs[2], 
                             mins[0], maxs[1], mins[2], 
                             mins[0], maxs[1], maxs[2],
                             maxs[0], mins[1], mins[2], 
                             maxs[0], mins[1], maxs[2], 
                             maxs[0], maxs[1], mins[2], 
                             maxs[0], maxs[1], maxs[2] };

  std::vector<IndexType> el { 3, 2, 4, 0, 
                              3, 4, 1, 0, 
                              3, 2, 6, 4,
                              3, 6, 7, 4, 
                              3, 5, 1, 4, 
                              3, 7, 5, 4 };
  // clang-format on

  points.swap(pt);
  elem.swap(el);
}

// 2D specialization for getBaryCoords(...)
template <>
Delaunay<2>::BaryCoordType Delaunay<2>::getBaryCoords(IndexType element_idx,
                                                      const PointType& query_pt) const
{
  const auto verts = m_mesh.ev_rel[element_idx];
  const ElementType tri(m_mesh.getVertexPoint(verts[0]),
                        m_mesh.getVertexPoint(verts[1]),
                        m_mesh.getVertexPoint(verts[2]));

  return tri.physToBarycentric(query_pt);
}

// 3D specialization for getBaryCoords(...)
template <>
Delaunay<3>::BaryCoordType Delaunay<3>::getBaryCoords(IndexType element_idx,
                                                      const PointType& query_pt) const
{
  const auto verts = m_mesh.ev_rel[element_idx];
  const ElementType tet(m_mesh.getVertexPoint(verts[0]),
                        m_mesh.getVertexPoint(verts[1]),
                        m_mesh.getVertexPoint(verts[2]),
                        m_mesh.getVertexPoint(verts[3]));

  return tet.physToBarycentric(query_pt);
}

// 2D specialization for isPointInSphere(...)
template <>
bool Delaunay<2>::InsertionHelper::isPointInSphere(const PointType& query_pt,
                                                   IndexType element_idx) const
{
  const auto verts = m_mesh.ev_rel[element_idx];
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  return primal::in_sphere(query_pt, p0, p1, p2, 0.);
}

// 3D specialization for isPointInSphere(...)
template <>
bool Delaunay<3>::InsertionHelper::isPointInSphere(const PointType& query_pt,
                                                   IndexType element_idx) const
{
  const auto verts = m_mesh.ev_rel[element_idx];
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  const PointType& p3 = m_mesh.getVertexPoint(verts[3]);
  return primal::in_sphere(query_pt, p0, p1, p2, p3, 0.);
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
