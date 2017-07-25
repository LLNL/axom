
/**
 * \file Delaunay.hpp
 *
 * \brief Construct delaunay triangulation by inserting points one by one.
 * A bounding box of the points needs to be defined first via initializeBoundary(...)
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
class Delaunay
{
public:
  using DataType = double;
  using Indextype = int;

#define TOPO_DIM 2
#define VERT_PER_ELEMENT 3
  using PointType = primal::Point<DataType, 2>;
  using IATriMeshType = slam::IAMesh<TOPO_DIM, 2, PointType>;
  using Triangle2D = axom::primal::Triangle<DataType, 2>;
  using BoundingBox = axom::primal::BoundingBox<DataType, 2>;

  using IndexListType = axom::slam::IndexListType;
  using IndexPairType = std::pair<IndexType, IndexType>;

  using CavityMapType = std::map<IndexListType, IndexListType>;
  using CavityMapPairType = std::pair<IndexListType, IndexListType>;

  using NewElemFacePairType = std::pair<IndexListType, IndexPairType>;
  using NewElemFaceMapType = std::map<IndexListType, IndexPairType>;

private:
  IATriMeshType m_mesh;
  BoundingBox m_bounding_box;
  bool m_has_boundary;
  int m_num_removed_elements_since_last_compact;

  CavityMapType m_cavity_boundary_map;

public:
  /**
   * \brief Default Constructor that creates the boundary square from -1,-1 to 1,1
   */
  Delaunay()
    : m_has_boundary(false)
    , m_num_removed_elements_since_last_compact(0)
  { }

  /**
   * \brief Defines the boundary of the triangulation.
   * \detail subsequent points added to the triangulation must not be outside of this boundary.
   */
  void initializeBoundary(const BoundingBox& bb)
  {
    m_bounding_box = bb;

    DataType xmin = bb.getMin()[0];
    DataType xmax = bb.getMax()[0];
    DataType ymin = bb.getMin()[1];
    DataType ymax = bb.getMax()[1];

    static const DataType point_arr[] =
      {xmin, ymin, xmin, ymax, xmax, ymin, xmax, ymax};
    std::vector<DataType> points(
      point_arr,
      point_arr + sizeof(point_arr) / sizeof(point_arr[0]));

    static const IndexType tri_arr[] = {0, 2, 1, 3, 1, 2};
    std::vector<IndexType> tri(tri_arr,
                               tri_arr + sizeof(tri_arr) / sizeof(tri_arr[0]));

    m_mesh = IATriMeshType(points, tri);

    m_has_boundary = true;
  }

  /**
   * \brief Add a new point to the m_mesh, which triggers re-triangulation
   * \brief This function will traverse the m_mesh to find the element that contains this point,
   * create the delaunay cavity, which takes out all the elements that contains the point it is sphere,
   * and fill in the cavity with a delaunay ball.
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

    //m_mesh.print_all();
    //SLIC_INFO("DT: Inserting point "<<new_pt[0] <<"," <<new_pt[1]);

    IndexType element_i = findContainingElement(new_pt);

    //SLIC_INFO("DT: Containing element: "<<element_i );

    std::vector<IndexType> elements_to_remove =
      findViolatingElements(new_pt, element_i);

    ///*
    std::cout << "DT: violating elements:";
    for(unsigned int i = 0; i < elements_to_remove.size(); i++)
      std::cout << elements_to_remove[i] << " ";
    std::cout << std::endl;
    //*/

    SLIC_ASSERT_MSG(elements_to_remove.size() > 0,
                    "Error: New point is not contained in the m_mesh");

    createCavity(elements_to_remove);

    //SLIC_INFO("DT: After removing elements");
    //m_mesh.print_all();

    //Add the new point
    IndexType new_pt_i = m_mesh.addVertex(new_pt);

    delaunayBall(new_pt_i);

    //SLIC_INFO("DT: After adding elements");
    //m_mesh.print_all();

    //Make sure ee_rel has only 4 invalid
    int invalid_nbr_count = 0;
    for(int i = 0; i < m_mesh.ee_rel.size(); i++)
    {
      if(m_mesh.ee_rel.isValidEntry(i))
      {
        for(int j = 0; j < VERT_PER_ELEMENT; j++)
        {
          invalid_nbr_count += m_mesh.ee_rel[i][j] < 0;
        }
      }
    }
    SLIC_ASSERT(invalid_nbr_count == 4);

    //call compact() if there are too many invalid points
    if(m_num_removed_elements_since_last_compact > 64 &&
       m_num_removed_elements_since_last_compact > (m_mesh.element_set.size() / 2))
    {
      m_mesh.compact();
      m_num_removed_elements_since_last_compact = 0;
    }

    if(!m_mesh.isValid(true))
    {
      SLIC_INFO("IA m_mesh is invalid after adding new point");
    }
  }

  /**
   * \brief Prints out the set, relation, and map detail of the m_mesh, for debug purpose.
   */
  void printMesh() { m_mesh.print_all(); }

  /**
   * \brief Write the m_mesh to a legacy VTK file
   * \param @filename the name of the file to write to, appending ".vtk" at the end of the filename
   * \detail This function uses mint to write the m_mesh to VTK format.
   */
  void writeToVTKFile(const std::string filename)
  {
    m_mesh.compact();

    axom::mint::UnstructuredMesh<MINT_TRIANGLE> mint_mesh(2);
    for(int i = 0; i < m_mesh.vertex_set.size(); i++)
    {
      if(!m_mesh.vertex_set.isValidEntry(i)) continue;
      mint_mesh.insertNode(m_mesh.vcoord_map[i][0], m_mesh.vcoord_map[i][1]);
    }
    for(int i = 0; i < m_mesh.ev_rel.size(); i++)
    {
      if(!m_mesh.ev_rel.isValidEntry(i)) continue;
      const int ptr[] = {m_mesh.ev_rel[i][0],
                         m_mesh.ev_rel[i][1],
                         m_mesh.ev_rel[i][2]};
      mint_mesh.insertCell(ptr, MINT_TRIANGLE, 3);
    }
    mint_mesh.toVtkFile(filename);
  }

  /**
    * \brief Removes the vertices that defines the boundary of the m_mesh, and the elements attached to them.
    * \detail After this function is called, no more points can be added to the m_mesh.
    */
  void removeBoundary()
  {
    if(m_has_boundary)
    {
      //remove the boundary box, which will be the first 4 points.
      m_mesh.removeVertex(0);
      m_mesh.removeVertex(1);
      m_mesh.removeVertex(2);
      m_mesh.removeVertex(3);

      m_has_boundary = false;
    }
  }

  /**
    * \brief Get the IA mesh data pointer
    */
  const IATriMeshType* getMeshData() const { return &m_mesh; }

private:
  Triangle2D getTriangle(IndexType element_idx)
  {
    std::vector<IndexType> verts = m_mesh.getVertexInElement(element_idx);

    Triangle2D tri(m_mesh.getVertexPoint(verts[0]),
                   m_mesh.getVertexPoint(verts[1]),
                   m_mesh.getVertexPoint(verts[2]));

    return tri;
  }

  /**
   * \brief Find the index of the element that contains the query point, or the element closest to the point.
   */
  IndexType findContainingElement(const PointType& query_pt)
  {
    SLIC_ASSERT(!m_mesh.isEmpty());

    //find the last valid element to use as starting element
    IndexType element_i = m_mesh.getValidElementIndex();

    SLIC_INFO("DT: find_containing_element " << element_i);
    SLIC_INFO("Query Pt " << query_pt);

    while(1)
    {
      SLIC_INFO("DT: step into element " << element_i);

      Triangle2D tri = getTriangle(element_i);

      axom::primal::Point3D bary_co = tri.barycentricCoords(query_pt);
      //SLIC_INFO("bary_co " << bary_co);

      //Find the most negative
      IndexType i = 0;
      if(bary_co[0] > bary_co[1]) i = 1;
      if(bary_co[i] > bary_co[2]) i = 2;

      if(bary_co[i] >= 0)
      {                    // bary_val bigger than zero -> inside triangle
        return element_i;  //return if inside or on triangle
      }

      std::vector<IndexType> zlist = m_mesh.getElementNeighbor(element_i);
      IndexType next_el = element_i = zlist[(i + 1) % 3];

      // Either there is a hole in the m_mesh, or the point is outside of the m_mesh.
      // Logically, this should never happen.
      SLIC_ASSERT(next_el >= 0);
    }
  }

  /**
   * \brief Find the list of element indices whose circle/sphere contains the query point.
   * \detail This function start from an element, and search through the neighboring elements,
   * returning a list of element indices whose circle/sphere defined by its vertices contains the
   * query point.
   * \param @query_pt the query point
   * \param @element_i the element to start the search at
   */
  std::vector<IndexType> findViolatingElements(const PointType& query_pt,
                                               IndexType element_i)
  {
    //SLIC_INFO("find_violating_elements from element " << q_element_i );

    std::vector<IndexType> ret;
    std::list<IndexType> element_list_to_check;
    element_list_to_check.push_back(element_i);

    std::set<IndexType> checked_elements;
    checked_elements.insert(element_i);
    checked_elements.insert(-1);

    //starting from element_i, which contains the point, try its neighbors, and the neighbor's neighbors
    //for point in circle test

    while(element_list_to_check.size() > 0)
    {
      element_i = element_list_to_check.front();
      element_list_to_check.pop_front();
      if(element_i < 0) continue;

      std::vector<IndexType> verts = m_mesh.getVertexInElement(element_i);
      const PointType p0 = m_mesh.getVertexPoint(verts[0]);
      const PointType p1 = m_mesh.getVertexPoint(verts[1]);
      const PointType p2 = m_mesh.getVertexPoint(verts[2]);
      bool is_in_circle = axom::primal::in_circle(query_pt, p0, p1, p2);

      if(is_in_circle)
      {
        ret.push_back(element_i);
        std::vector<IndexType> nbr_elem_list =
          m_mesh.getElementNeighbor(element_i);
        for(int i = 0; i < (int)nbr_elem_list.size(); i++)
        {
          IndexType nbr_element = nbr_elem_list[i];
          std::set<IndexType>::iterator it = checked_elements.find(nbr_element);
          if(it == checked_elements.end())
          {
            element_list_to_check.push_back(nbr_element);
            checked_elements.insert(nbr_element);
          }
        }
      }
    }
    return ret;
  }

  /**
    * \brief Remove the elements in the delaunay cavity
    */
  void createCavity(const std::vector<IndexType>& elements_to_remove)
  {
    //m_cavity_boundary_map: <vlist, <vlist, element_idx> >
    // m_cavity_boundary_map is a map from sorted face vertices to pair of face vertices and element idx.
    // It keeps track of the edges that forms the cavity, which will be edges of the new element.
    m_cavity_boundary_map.clear();

    // Take out the elements one by one, keeping track of the cavity edges
    for(unsigned int i = 0; i < elements_to_remove.size(); i++)
    {
      IndexType element_i = elements_to_remove[i];

      //For each face of the element, insert its face into the map (or remove if duplicate)
      for(unsigned int face_i = 0; face_i < VERT_PER_ELEMENT; face_i++)
      {
        IndexListType face_vlist = m_mesh.getElementFace(element_i, face_i);

        IndexListType face_vlist_sorted(
          face_vlist);  //sort the list to make the key for the map
        std::sort(face_vlist_sorted.begin(), face_vlist_sorted.end());

        std::pair<CavityMapType::iterator, bool> ret =
          m_cavity_boundary_map.insert(
            CavityMapPairType(face_vlist_sorted, face_vlist));
        if(!ret.second)
        {
          m_cavity_boundary_map.erase(ret.first);
        }
      }

      m_mesh.removeElement(element_i);
      m_num_removed_elements_since_last_compact++;
    }
  }

  /**
    * \brief Fill in the delaunay cavity
    */
  void delaunayBall(IndexType new_pt_i)
  {
    //new_element_face_map: <vlist, <element_idx,face_idx> >
    // new_element_face_map is a map from sorted face vertices to pair of element and its face idx.
    // It keep tracks of the new elements added, to correctly construct element->element relation
    NewElemFaceMapType new_element_face_map;

    //Add new triangles from the cavity edges
    for(CavityMapType::iterator cav_edge_iter = m_cavity_boundary_map.begin();
        cav_edge_iter != m_cavity_boundary_map.end();
        cav_edge_iter++)
    {
      IndexType n0 = cav_edge_iter->second[0];
      IndexType n1 = cav_edge_iter->second[1];
      IndexType new_el = m_mesh.addElement(n0, n1, new_pt_i);
      //SLIC_INFO("new element " << new_el <<": "<< n0 <<" "<< n1<< " " <<new_pt_i);

      // Check the new element's face neighbor, because it can be wrong sometimes
      // due to the m_mesh being temporarily non-manifold.
      for(int fi = 0; fi < VERT_PER_ELEMENT; fi++)
      {
        IndexListType vlist_sorted = m_mesh.getElementFace(new_el, fi);
        std::sort(vlist_sorted.begin(), vlist_sorted.end());

        NewElemFaceMapType::iterator iter =
          new_element_face_map.find(vlist_sorted);

        if(iter == new_element_face_map.end())
        {
          new_element_face_map.insert(
            NewElemFacePairType(vlist_sorted, IndexPairType(new_el, fi)));
        }
        else
        {
          IndexType nbr_el = iter->second.first;
          IndexType nbr_fi = iter->second.second;

          if(m_mesh.ee_rel[new_el][fi] < 0)
          {
            //SLIC_INFO("found neighbor " << iter->second.first <<"-"<<nbr_fi);

            m_mesh.ee_rel.modify(new_el, fi, nbr_el);
            m_mesh.ee_rel.modify(nbr_el, nbr_fi, new_el);
          }

          new_element_face_map.erase(iter);
        }
      }
    }
  }
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
