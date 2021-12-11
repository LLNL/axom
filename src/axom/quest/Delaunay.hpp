
#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

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
  using IAMeshType = slam::IAMesh<DIM, DIM, PointType>;

  using Indextype = typename IAMeshType::IndexType;
  using IndexArray = typename IAMeshType::IndexArray;

  using BaryCoordType = primal::Point<DataType, DIM + 1>;
  using Point2DType = primal::Point<DataType, 2>;
  using Point3DType = primal::Point<DataType, 3>;
  using Point4DType = primal::Point<DataType, 4>;
  using Triangle2D = axom::primal::Triangle<DataType, 2>;
  using Tetrahedron3D = axom::primal::Tetrahedron<DataType, 3>;
  using BoundingBox = axom::primal::BoundingBox<DataType, DIM>;

  using IndexPairType = std::pair<IndexType, IndexType>;

  static constexpr int VERT_PER_ELEMENT = DIM + 1;
  static constexpr IndexType INVALID_ID = -1;

private:
  using ModularFaceIndex = axom::slam::ModularInt<
    axom::slam::policies::CompileTimeSize<IndexType, VERT_PER_ELEMENT>>;

  /// Helper struct for storing the facet indices of each simplex
  template <int TOP_DIM>
  struct ElementFacePair
  {
    IndexType element_idx;
    IndexType face_vidx[TOP_DIM];
    ElementFacePair(IndexType el, IndexType* vlist)
    {
      element_idx = el;
      for(int i = 0; i < TOP_DIM; i++)
      {
        face_vidx[i] = vlist[i];
      }
    }
  };

private:
  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;

  std::vector<ElementFacePair<DIM>> cavity_face_list;
  IndexArray cavity_element_list;
  IndexArray new_elements;
  std::set<IndexType> checked_element_set;

public:
  /**
   * \brief Default Constructor. User need to call initializeBoundary(BoundingBox) before adding points.
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
   * \brief Adds a new point to the mesh and locally re-triangulates
   * the mesh to ensure that it stays Delaunay.
   *
   * This function will traverse the m_mesh to find the element that contains
   * this point, create the Delaunay cavity, which takes out all the elements
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

    IndexType element_i = findContainingElement(new_pt);

    if(element_i == INVALID_ID)
    {
      SLIC_WARNING(axom::fmt::format(
        "Could not insert point {} into Delaunay triangulation: Element "
        "containing that point was not found"));
      return;
    }

    findCavityElements(new_pt, element_i);

    createCavity();

    //Add the new point
    IndexType new_pt_i = m_mesh.addVertex(new_pt);

    delaunayBall(new_pt_i);

    m_mesh.fixVertexNeighborhood(new_pt_i, new_elements);

    //call compact() if there are too many invalid points
    //This auto-compacting feature is hard coded. It may be good to let user have control
    // of this option in the future.
    if(m_num_removed_elements_since_last_compact > 64 &&
       m_num_removed_elements_since_last_compact > (m_mesh.element_set.size() / 2))
    {
      m_mesh.compact();
      m_num_removed_elements_since_last_compact = 0;
    }

    SLIC_WARNING_IF(!m_mesh.isValid(),
                    "IA m_mesh is invalid after adding new point");
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
    m_mesh.compact();
    const auto CELL_TYPE = DIM == 2 ? mint::TRIANGLE : mint::TET;

    using UMesh = mint::UnstructuredMesh<mint::SINGLE_SHAPE>;
    UMesh mint_mesh(DIM, CELL_TYPE);

    for(int i = 0; i < m_mesh.vertex_set.size(); i++)
    {
      mint_mesh.appendNodes(m_mesh.getVertexPoint(i).data(), 1);
    }

    for(int i = 0; i < m_mesh.ev_rel.size(); i++)
    {
      const auto* ptr = &(m_mesh.ev_rel[i][0]);
      mint_mesh.appendCell(ptr, CELL_TYPE);
    }
    mint::write_vtk(&mint_mesh, filename);
  }

  /**
   * \brief Removes the vertices that defines the boundary of the mesh,
   * and the elements attached to them.
   *
   * \details After this function is called,
   * no more points can be added to the m_mesh.
   */
  void removeBoundary()
  {
    if(m_has_boundary)
    {
      //remove the boundary box, which will be the first 4 points for triangles, first 8 for tetrahedron
      const int num_boundary_pts = DIM == 2 ? 4 : 8;

      //Collect a list of elements to remove first, because
      //the list may be incomplete if generated during the removal.
      IndexArray elements_to_remove;
      for(int i = 0; i < num_boundary_pts; i++)
      {
        IndexArray elist = m_mesh.getElementsWithVertex(i);
        elements_to_remove.insert(elements_to_remove.end(),
                                  elist.begin(),
                                  elist.end());
      }
      for(unsigned int i = 0; i < elements_to_remove.size(); i++)
      {
        if(m_mesh.isValidElementEntry(elements_to_remove[i]))
          m_mesh.removeElement(elements_to_remove[i]);
      }

      for(int i = 0; i < num_boundary_pts; i++)
      {
        m_mesh.removeVertex(i);
      }

      m_has_boundary = false;
    }
  }

  /// \brief Get the IA mesh data pointer
  const IAMeshType* getMeshData() const { return &m_mesh; }

private:
  /// \brief Find the index of the element that contains the query point, or the element closest to the point.
  IndexType findContainingElement(const PointType& query_pt)
  {
    if(m_mesh.isEmpty())
    {
      SLIC_ERROR(
        "Attempting to insert point into empty Delaunay triangulation."
        "Delaunay::initializeBoundary() needs to be called first");
      return INVALID_ID;
    }
    if(!m_bounding_box.contains(query_pt))
    {
      SLIC_WARNING(
        "Attempting to locate element at location outside valid domain");
      return INVALID_ID;
    }

    //find the last valid element to use as starting element
    IndexType element_i = m_mesh.getValidElementIndex();

    while(1)
    {
      BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);

      //Find the index of the most negative barycentric coord
      ModularFaceIndex idx(bary_coord.array().argMin());

      if(bary_coord[idx] >= 0)  // inside if smallest bary coord positive
      {
        return element_i;
      }

      // else, move to that neighbor
      IndexArray zlist = m_mesh.getElementNeighbors(element_i);
      element_i = zlist[idx + 1];

      // Either there is a hole in the m_mesh, or the point is outside of the m_mesh.
      // Logically, this should never happen.
      SLIC_ASSERT(m_mesh.isValidElementEntry(element_i));
    }
  }

  /**
   * \brief Helper function returns true if the query point is in the sphere formed by the element vertices
   */
  bool isPointInSphere(const PointType& query_pt, IndexType element_idx);

  /**
   * \brief recursive function to find cavity elements given a point to be added
   * \details Check if the point is in the circle/sphere of the element, if so,
   * recursively call the neighboring elements.
   */
  bool findCavityElementsRec(const PointType& query_pt, IndexType element_idx)
  {
    const bool is_in_circle = isPointInSphere(query_pt, element_idx);

    if(is_in_circle)
    {
      //add to cavity elements
      cavity_element_list.push_back(element_idx);

      //check for each faces
      IndexArray nbr_elements = m_mesh.getElementNeighbors(element_idx);
      SLIC_ASSERT(nbr_elements.size() == VERT_PER_ELEMENT);

      for(int face_i = 0; face_i < (int)VERT_PER_ELEMENT; face_i++)
      {
        IndexType nbr_elem = nbr_elements[face_i];

        //The latter case is a checked element that is not a cavity element
        if(!m_mesh.isValidElementEntry(nbr_elem) ||
           (checked_element_set.insert(nbr_elem).second
              ? findCavityElementsRec(query_pt, nbr_elem)
              : !axom::slam::is_subset(nbr_elem, cavity_element_list)))
        {
          IndexArray vlist = m_mesh.getElementFace(element_idx, face_i);

          //For tetrahedron, if the element face is odd, reverse vertex order
          if(DIM == 3 && face_i % 2 == 1)
          {
            axom::utilities::swap(vlist[1], vlist[2]);
          }

          cavity_face_list.push_back(ElementFacePair<DIM>(element_idx, &vlist[0]));
        }
      }
      return false;
    }
    else
    {
      //add to cavity faces
      return true;
    }
  }

  /**
   * \brief Find the list of element indices whose sphere contains query point
   *
   * \details This function start from an element, and search through the
   * neighboring elements, returning a list of element indices
   * whose circle/sphere defined by its vertices contains the query point.
   * \param query_pt the query point
   * \param element_i the element to start the search at
   */
  void findCavityElements(const PointType& query_pt, IndexType element_i)
  {
    cavity_element_list.clear();
    cavity_face_list.clear();
    checked_element_set.clear();
    checked_element_set.insert(element_i);
    findCavityElementsRec(query_pt, element_i);

    SLIC_ASSERT_MSG(cavity_element_list.size() > 0,
                    "Error: New point is not contained in the mesh");
    SLIC_ASSERT(cavity_face_list.size() > 0);
  }

  /**
    * \brief Remove the elements in the delaunay cavity
    */
  void createCavity()
  {
    for(unsigned int i = 0; i < cavity_element_list.size(); i++)
    {
      m_mesh.removeElement(cavity_element_list[i]);
      ++m_num_removed_elements_since_last_compact;
    }
  }

  /// \brief Fill in the delaunay cavity
  void delaunayBall(IndexType new_pt_i)
  {
    new_elements.clear();

    for(unsigned int i = 0; i < cavity_face_list.size(); i++)
    {
      IndexType vlist[VERT_PER_ELEMENT];
      for(int d = 0; d < VERT_PER_ELEMENT - 1; ++d)
      {
        vlist[d] = cavity_face_list[i].face_vidx[d];
      }
      vlist[VERT_PER_ELEMENT - 1] = new_pt_i;

      IndexType new_el = m_mesh.addElement(vlist);
      new_elements.push_back(new_el);
    }
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
  BaryCoordType getBaryCoords(IndexType element_idx, const PointType& q_pt);

};  //END class Delaunay

//********************************************************************************
// Below are 2D and 3D specializations for methods in the Delaunay class
//********************************************************************************

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
                                                      const PointType& query_pt)
{
  IndexArray verts = m_mesh.getVerticesInElement(element_idx);

  Triangle2D tri(m_mesh.getVertexPoint(verts[0]),
                 m_mesh.getVertexPoint(verts[1]),
                 m_mesh.getVertexPoint(verts[2]));

  BaryCoordType bary_co = tri.physToBarycentric(query_pt);

  return bary_co;
}

// 3D specialization for getBaryCoords(...)
template <>
Delaunay<3>::BaryCoordType Delaunay<3>::getBaryCoords(IndexType element_idx,
                                                      const PointType& query_pt)
{
  IndexArray verts = m_mesh.getVerticesInElement(element_idx);

  Tetrahedron3D tet(m_mesh.getVertexPoint(verts[0]),
                    m_mesh.getVertexPoint(verts[1]),
                    m_mesh.getVertexPoint(verts[2]),
                    m_mesh.getVertexPoint(verts[3]));

  BaryCoordType bary_co = tet.physToBarycentric(query_pt);

  return bary_co;
}

// 2D specialization for isPointInSphere(...)
template <>
bool Delaunay<2>::isPointInSphere(const PointType& query_pt, IndexType element_idx)
{
  IndexArray verts = m_mesh.getVerticesInElement(element_idx);
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  return axom::primal::in_sphere(query_pt, p0, p1, p2);
}

// 3D specialization for isPointInSphere(...)
template <>
bool Delaunay<3>::isPointInSphere(const PointType& query_pt, IndexType element_idx)
{
  IndexArray verts = m_mesh.getVerticesInElement(element_idx);
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  const PointType& p3 = m_mesh.getVertexPoint(verts[3]);
  return axom::primal::in_sphere(query_pt, p0, p1, p2, p3);
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
