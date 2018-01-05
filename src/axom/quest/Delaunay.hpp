
/**
 * \file Delaunay.hpp
 *
 * \brief Construct delaunay triangulation by inserting points one by one.
 * A bounding box of the points needs to be defined first via initializeBoundary(...)
 *
 * A bounding box of the points needs to be defined first via
 * the initializeBoundary() function
 */

#ifndef QUEST_DELAUNAY_H_
#define QUEST_DELAUNAY_H_

#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"  //for writing out to VTK file

#include <list>
#include <vector>
#include <set>
#include <stdlib.h>

namespace axom
{
namespace quest
{
template <unsigned int DIMENSION = 2>
class Delaunay
{
public:
  AXOM_STATIC_ASSERT_MSG(DIMENSION == 2 || DIMENSION == 3,
                         "The template parameter DIMENSION can only be 2 or 3. "
                         "Other dimensions are not supported");

  enum
  {
    VERT_PER_ELEMENT = DIMENSION + 1
  };

  using DataType = double;
  using Indextype = int;

  using PointType = primal::Point<DataType, DIMENSION>;
  using BaryCoordType = primal::Point<DataType, DIMENSION + 1>;
  using Point2DType = primal::Point<DataType, 2>;
  using Point3DType = primal::Point<DataType, 3>;
  using Point4DType = primal::Point<DataType, 4>;
  using IAMeshType = slam::IAMesh<DIMENSION, DIMENSION, PointType>;
  using Triangle2D = axom::primal::Triangle<DataType, 2>;
  using Tetrahedron3D = axom::primal::Tetrahedron<DataType, 3>;
  using BoundingBox = axom::primal::BoundingBox<DataType, DIMENSION>;

  using IndexListType = typename IAMeshType::IndexListType;
  using IndexPairType = std::pair<IndexType, IndexType>;

private:
  IAMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;

  template <unsigned int TOP_DIM>
  struct ElementFacePair
  {
    int element_idx;
    int face_vidx[TOP_DIM];
    ElementFacePair(int el, int* vlist)
    {
      element_idx = el;
      for(int i = 0; i < (int)TOP_DIM; i++) face_vidx[i] = vlist[i];
    }
  };

  std::vector<ElementFacePair<DIMENSION>> cavity_face_list;
  std::vector<IndexType> cavity_element_list;
  std::vector<IndexType> new_elements;
  std::set<IndexType> checked_element_set;

public:
  /**
   * \brief Default Constructor. User need to call initializeBoundary(BoundingBox) before adding points.
   */
  Delaunay()
    : m_has_boundary(false)
    , m_num_removed_elements_since_last_compact(0)
  {
    SLIC_ASSERT(DIMENSION == 2 || DIMENSION == 3);
  }

  /**
   * \brief Defines the boundary of the triangulation.
   * \details subsequent points added to the triangulation must not be outside of this boundary.
   */
  void initializeBoundary(const BoundingBox& bb)
  {
    std::vector<DataType> points;
    std::vector<IndexType> elem;

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

    SLIC_WARNING_IF(!m_mesh.isValid(true),
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
   * \note The suffix ".vtk" will be appended to the provided  filename
   * \details This function uses mint to write the m_mesh to VTK format.
   */
  void writeToVTKFile(const std::string filename)
  {
    m_mesh.compact();
    const auto CELL_TYPE =
      DIMENSION == 2 ? axom::mint::TRIANGLE : axom::mint::TET;

    axom::mint::UnstructuredMesh<mint::SINGLE_SHAPE> mint_mesh(DIMENSION,
                                                               CELL_TYPE);

    for(int i = 0; i < m_mesh.vertex_set.size(); i++)
    {
      mint_mesh.appendNodes(m_mesh.getVertexPoint(i).data());
    }

    for(int i = 0; i < m_mesh.ev_rel.size(); i++)
    {
      const int* ptr = &(m_mesh.ev_rel[i][0]);
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
      unsigned int num_boundary_pts = DIMENSION == 2 ? 4 : 8;

      //Collect a list of elements to remove first, because
      //the list may be incomplete if generated during the removal.
      IndexListType elements_to_remove;
      for(unsigned int i = 0; i < num_boundary_pts; i++)
      {
        IndexListType elist = m_mesh.getElementsWithVertex(i);
        elements_to_remove.insert(elements_to_remove.end(),
                                  elist.begin(),
                                  elist.end());
      }
      for(unsigned int i = 0; i < elements_to_remove.size(); i++)
      {
        if(m_mesh.isValidElementEntry(elements_to_remove[i]))
          m_mesh.removeElement(elements_to_remove[i]);
      }

      for(unsigned int i = 0; i < num_boundary_pts; i++) m_mesh.removeVertex(i);

      m_has_boundary = false;
    }
  }

  /**
    * \brief Get the IA mesh data pointer
    */
  const IAMeshType* getMeshData() const { return &m_mesh; }

private:
  /**
   * \brief Find the index of the element that contains the query point, or the element closest to the point.
   */
  IndexType findContainingElement(const PointType& query_pt)
  {
    SLIC_ASSERT(!m_mesh.isEmpty());

    //find the last valid element to use as starting element
    IndexType element_i = m_mesh.getValidElementIndex();

    while(1)
    {
      BaryCoordType bary_coord = getBaryCoords(element_i, query_pt);

      //Find the most negative
      IndexType i = 0;
      for(unsigned int d = 1; d < VERT_PER_ELEMENT; d++)
      {
        if(bary_coord[i] > bary_coord[d]) i = d;
      }

      if(bary_coord[i] >= 0)
      {                    // bary_val bigger than zero -> inside triangle
        return element_i;  //return if inside or on triangle
      }

      std::vector<IndexType> zlist = m_mesh.getElementNeighbors(element_i);
      element_i = zlist[(i + 1) % VERT_PER_ELEMENT];

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
    bool is_in_circle = isPointInSphere(query_pt, element_idx);

    if(is_in_circle)
    {
      //add to cavity elements
      cavity_element_list.push_back(element_idx);

      //check for each faces
      IndexListType nbr_elements = m_mesh.getElementNeighbors(element_idx);
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
          IndexListType vlist = m_mesh.getElementFace(element_idx, face_i);

          //For tetrahedron, if the element face is odd, reverse vertex order
          if(DIMENSION == 3 && face_i % 2 == 1)
          {
            IndexType tmp = vlist[1];
            vlist[1] = vlist[2];
            vlist[2] = tmp;
          }

          cavity_face_list.push_back(
            ElementFacePair<DIMENSION>(element_idx, &vlist[0]));
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

    return;
  }

  /**
    * \brief Remove the elements in the delaunay cavity
    */
  void createCavity()
  {
    for(unsigned int i = 0; i < cavity_element_list.size(); i++)
    {
      m_mesh.removeElement(cavity_element_list[i]);
      m_num_removed_elements_since_last_compact++;
    }

    return;
  }

  /**
    * \brief Fill in the delaunay cavity
    */
  void delaunayBall(IndexType new_pt_i)
  {
    new_elements.clear();

    for(unsigned int i = 0; i < cavity_face_list.size(); i++)
    {
      IndexType vlist[VERT_PER_ELEMENT];
      for(unsigned int d = 0; d < VERT_PER_ELEMENT - 1; d++)
      {
        vlist[d] = cavity_face_list[i].face_vidx[d];
      }
      vlist[VERT_PER_ELEMENT - 1] = new_pt_i;

      IndexType new_el = m_mesh.addElement(vlist);
      new_elements.push_back(new_el);
    }

    return;
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

// this is the 2D specialization for generateInitialMesh(...)
template <>
void Delaunay<2>::generateInitialMesh(std::vector<DataType>& points,
                                      std::vector<IndexType>& elem,
                                      const BoundingBox& bb)
{
  //Set up the initial IA mesh of 2 triangles forming a rectangle

  DataType xmin = bb.getMin()[0];
  DataType xmax = bb.getMax()[0];
  DataType ymin = bb.getMin()[1];
  DataType ymax = bb.getMax()[1];

  static const DataType point_arr[] = {xmin, ymin, xmin, ymax, xmax, ymin, xmax, ymax};
  static const IndexType tri_arr[] = {0, 2, 1, 3, 1, 2};
  points.assign(point_arr, point_arr + sizeof(point_arr) / sizeof(point_arr[0]));
  elem.assign(tri_arr, tri_arr + sizeof(tri_arr) / sizeof(tri_arr[0]));
}

// this is the 3D specialization for generateInitialMesh(...)
template <>
void Delaunay<3>::generateInitialMesh(std::vector<DataType>& points,
                                      std::vector<IndexType>& elem,
                                      const BoundingBox& bb)
{
  //Set up the initial IA mesh of 6 tetrahedrons forming a cube

  DataType xmin = bb.getMin()[0];
  DataType xmax = bb.getMax()[0];
  DataType ymin = bb.getMin()[1];
  DataType ymax = bb.getMax()[1];
  DataType zmin = bb.getMin()[2];
  DataType zmax = bb.getMax()[2];

  static const DataType point_arr[] = {
    xmin, ymin, zmin, xmin, ymin, zmax, xmin, ymax, zmin, xmin, ymax, zmax,
    xmax, ymin, zmin, xmax, ymin, zmax, xmax, ymax, zmin, xmax, ymax, zmax};

  static const IndexType tet_arr[] = {3, 2, 4, 0, 3, 4, 1, 0, 3, 2, 6, 4,
                                      3, 6, 7, 4, 3, 5, 1, 4, 3, 7, 5, 4};

  points.assign(point_arr, point_arr + sizeof(point_arr) / sizeof(point_arr[0]));
  elem.assign(tet_arr, tet_arr + sizeof(tet_arr) / sizeof(tet_arr[0]));
}

// this is the 2D specialization for getBaryCoords(...)
template <>
Delaunay<2>::BaryCoordType Delaunay<2>::getBaryCoords(IndexType element_idx,
                                                      const PointType& query_pt)
{
  IndexListType verts = m_mesh.getVerticesInElement(element_idx);

  Triangle2D tri(m_mesh.getVertexPoint(verts[0]),
                 m_mesh.getVertexPoint(verts[1]),
                 m_mesh.getVertexPoint(verts[2]));

  BaryCoordType bary_co = tri.physToBarycentric(query_pt);

  return bary_co;
}

// this is the 3D specialization for getBaryCoords(...)
template <>
Delaunay<3>::BaryCoordType Delaunay<3>::getBaryCoords(IndexType element_idx,
                                                      const PointType& query_pt)
{
  IndexListType verts = m_mesh.getVerticesInElement(element_idx);

  Tetrahedron3D tet(m_mesh.getVertexPoint(verts[0]),
                    m_mesh.getVertexPoint(verts[1]),
                    m_mesh.getVertexPoint(verts[2]),
                    m_mesh.getVertexPoint(verts[3]));

  BaryCoordType bary_co = tet.physToBarycentric(query_pt);

  return bary_co;
}

// this is the 2D specialization for isPointInSphere(...)
template <>
bool Delaunay<2>::isPointInSphere(const PointType& query_pt, IndexType element_idx)
{
  IndexListType verts = m_mesh.getVerticesInElement(element_idx);
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  return axom::primal::in_sphere(query_pt, p0, p1, p2);
}

// this is the 3D specialization for isPointInSphere(...)
template <>
bool Delaunay<3>::isPointInSphere(const PointType& query_pt, IndexType element_idx)
{
  IndexListType verts = m_mesh.getVerticesInElement(element_idx);
  const PointType& p0 = m_mesh.getVertexPoint(verts[0]);
  const PointType& p1 = m_mesh.getVertexPoint(verts[1]);
  const PointType& p2 = m_mesh.getVertexPoint(verts[2]);
  const PointType& p3 = m_mesh.getVertexPoint(verts[3]);
  return axom::primal::in_sphere(query_pt, p0, p1, p2, p3);
}

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
