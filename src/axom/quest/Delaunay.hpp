
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
#define TOPO_DIM 2
#define VERT_PER_ELEMENT TOPO_DIM + 1
#define SPATIAL_DIM 2

  using DataType = double;
  using Indextype = int;

  using PointType = primal::Point<DataType, SPATIAL_DIM>;
  using IATriMeshType = slam::IAMesh<TOPO_DIM, SPATIAL_DIM, PointType>;
  using Triangle2D = axom::primal::Triangle<DataType, SPATIAL_DIM>;
  using BoundingBox = axom::primal::BoundingBox<DataType, SPATIAL_DIM>;

  using IndexListType = IATriMeshType::IndexListType;
  using IndexPairType = std::pair<IndexType, IndexType>;

private:
  IATriMeshType m_mesh;
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

  std::vector<ElementFacePair<TOPO_DIM>> cavity_face_list;
  std::vector<IndexType> cavity_element_list;
  std::vector<IndexType> new_elements;
  std::set<IndexType> checked_element_set;

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

    findCavityElements(new_pt, element_i);

    /*
    std::cout<<"DT: violating elements:" ;
    for(unsigned int i=0; i<elements_to_remove.size(); i++)
      std::cout<<elements_to_remove[i] <<" ";
    std::cout<<std::endl;
    */

    //SLIC_ASSERT_MSG(elements_to_remove.size() > 0, "Error: New point is not contained in the m_mesh");

    createCavity();

    //SLIC_INFO("DT: After removing elements");
    //m_mesh.print_all();

    //Add the new point
    IndexType new_pt_i = m_mesh.addVertex(new_pt);

    delaunayBall(new_pt_i);

    m_mesh.fixVertexNeighborhood(new_pt_i, new_elements);

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

      const int* ptr = &(m_mesh.ev_rel[i][0]);

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
    IndexListType verts = m_mesh.getVertexInElement(element_idx);

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

    //SLIC_INFO("DT: find_containing_element " << element_i );
    //SLIC_INFO("Query Pt " << query_pt);

    while(1)
    {
      //SLIC_INFO("DT: step into element "<<element_i);

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
   * \brief recursive function to find cavity elements given a point to be added
   * \detail Check if the point is in the circle/sphere of the element, if so, call
   * recursively on the neighboring elements.
   */
  bool findCavityElementsRec(const PointType& query_pt, IndexType element_idx)
  {
    IndexListType verts = m_mesh.getVerticesInElement(element_idx);
    const PointType p0 = m_mesh.getVertexPoint(verts[0]);
    const PointType p1 = m_mesh.getVertexPoint(verts[1]);
    const PointType p2 = m_mesh.getVertexPoint(verts[2]);
    bool is_in_circle = axom::primal::in_circle(query_pt, p0, p1, p2);

    //SLIC_INFO(element_idx << " InCircle:" << is_in_circle);

    if(is_in_circle)
    {
      //add to cavity elements
      cavity_element_list.push_back(element_idx);

      //check for each faces
      IndexListType nbr_elements = m_mesh.getElementNeighbors(element_idx);
      SLIC_ASSERT(nbr_elements.size() == TOPO_DIM + 1);

      for(int face_i = 0; face_i < VERT_PER_ELEMENT; face_i++)
      {
        IndexType nbr_elem = nbr_elements[face_i];

        //SLIC_INFO("nbr "<< face_i <<" "<<nbr_elem);

        if(nbr_elem < 0 ||
           (checked_element_set.insert(nbr_elem).second
              ? findCavityElementsRec(query_pt, nbr_elem)
              : !axom::slam::is_subset(
                  nbr_elem,
                  cavity_element_list)  //checked but not a cavity element
            ))
        {
          IndexListType vlist = m_mesh.getElementFace(element_idx, face_i);
          cavity_face_list.push_back(
            ElementFacePair<TOPO_DIM>(element_idx, &vlist[0]));
          //SLIC_INFO("Added "<<vlist[0]<< " "<<vlist[1]);
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
   * \brief Find the list of element indices whose circle/sphere contains the query point.
   * \detail This function start from an element, and search through the neighboring elements,
   * returning a list of element indices whose circle/sphere defined by its vertices contains the
   * query point.
   * \param @query_pt the query point
   * \param @element_i the element to start the search at
   */
  void findCavityElements(const PointType& query_pt, IndexType element_i)
  {
    //SLIC_INFO("find_violating_elements from element " << element_i );

    cavity_element_list.clear();
    cavity_face_list.clear();
    checked_element_set.clear();
    checked_element_set.insert(element_i);
    findCavityElementsRec(query_pt, element_i);

    /*
    std::cout<<"violating elements: " ;
    for(unsigned int i=0; i<cavity_element_list.size(); i++)
    {
      std::cout<< cavity_element_list[i] <<" ";
    }
    std::cout<< std::endl;

    std::cout<<"adding elements: " ;
    for(unsigned int i=0; i<cavity_face_list.size(); i++)
    {
      std::cout<< cavity_face_list[i].face_vidx[0] << "-" << cavity_face_list[i].face_vidx[1] <<", ";
    }
    std::cout<< std::endl;
    //*/

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
      IndexType n0 = cavity_face_list[i].face_vidx[0];
      IndexType n1 = cavity_face_list[i].face_vidx[1];
      IndexType new_el = m_mesh.addElement(n0, n1, new_pt_i);
      new_elements.push_back(new_el);
    }

    return;
  }
};

}  // end namespace quest
}  // end namespace axom

#endif  //  QUEST_DELAUNAY_H_
