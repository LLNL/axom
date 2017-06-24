// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*
 * \file TopologicalMesh.cpp
 *
 * \brief Example of using SLAM to build topological data structure
 *
 */

#include "axom/core.hpp"
#include "axom/slam.hpp"
#include "axom/primal.hpp"

#include <list>
#include <set>
#include <stdlib.h>

namespace slamTopologicalMesh
{
/// types for sets
using NodeSet = axom::slam::PositionSet<>;
using ZoneSet = axom::slam::PositionSet<>;
using PositionType = ZoneSet::PositionType;
using IndexType = ZoneSet::ElementType;
using DataType = double;

struct BasicTriMeshData
{
  //This is a regular cube
  std::vector<double> points = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                                0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 1.0,
                                0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0};
  std::vector<int> tri = {0, 6, 4, 0, 2, 6, 0, 3, 2, 0, 1, 3, 2, 7, 6, 2, 3, 7,
                          4, 6, 7, 4, 7, 5, 0, 4, 5, 0, 5, 1, 1, 5, 7, 1, 7, 3};
};

struct Basic2DTriMeshData
{
  //This is a square
  std::vector<double> points = {0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0};
  std::vector<int> tri = {0, 1, 2, 1, 3, 2};
};

struct BasicTetMeshData
{
  //cube divided into tets

  std::vector<double> points = {-1.0, -1.0, -1.0, -1.0, -1.0, 1.0,  -1.0, 1.0,
                                -1.0, -1.0, 1.0,  1.0,  1.0,  -1.0, -1.0, 1.0,
                                -1.0, 1.0,  1.0,  1.0,  -1.0, 1.0,  1.0,  1.0};

  std::vector<int> tet = {3, 2, 4, 0, 3, 1, 4, 0, 3, 6, 2, 4,
                          3, 6, 7, 4, 3, 5, 1, 4, 3, 5, 7, 4};
};

struct Flat2DTriMeshData
{
  //This is a square
  std::vector<double> points =
    {-1, -1, -1, 1, 1, -1, 1, 1, -2, -2, -2, 2, 2, -2, 2, 2};
  std::vector<int> tri = {

    1, 3, 5, 4, 0, 5, 0, 1, 5, 0, 2, 1, 0, 4, 2,
    3, 7, 5, 3, 2, 6, 3, 1, 2, 3, 6, 7, 2, 4, 6};
};

typedef axom::primal::Point2D Point2D;
typedef axom::primal::Point3D Point3D;
typedef axom::primal::Vector2D Vector2D;
typedef axom::primal::Triangle<double, 2> Triangle2D;
//typedef axom::slam::Map< Point2D >            XYZField;
/*
	void calculate_normal(axom::slam::IAMesh<3> & tri_mesh){

		//calculate normal of each zone as a map
		XYZField zone_to_normal_map = XYZField(&tri_mesh.zone_set);
		for(IndexType zidx = 0; zidx < tri_mesh.zone_set.size(); zidx++)
		{

			Point p0 = tri_mesh.dcoord_map[ tri_mesh.dzn_rel[zidx][0] ];
			Point p1 = tri_mesh.dcoord_map[ tri_mesh.dzn_rel[zidx][1] ];
			Point p2 = tri_mesh.dcoord_map[ tri_mesh.dzn_rel[zidx][2] ];
			Point p01 = p1 - p0;
			Point p02 = p2 - p0;
			//Do we assume the triangle is not degenerate?
			Point pnormal = cross( p01, p02);
			zone_to_normal_map[zidx] = pnormal;
			//SLIC_INFO("normal "<< pnormal.m_x <<" " << pnormal.m_y << " " <<pnormal.m_z);
		}

		//calculate normal of each node as a map
		XYZField node_to_normal_map = XYZField(&tri_mesh.zone_set);
		for(IndexType nidx = 0; nidx < tri_mesh.node_set.size(); nidx++)
		{
			Point node_normal(0,0,0);
			std::vector<IndexType> zone_subset = tri_mesh.getZoneWithNode(nidx);
			for(unsigned int zidx = 0; zidx < zone_subset.size(); zidx++)
			{
				Point& zone_normal = zone_to_normal_map[zone_subset[zidx]];
				node_normal += zone_normal;
			}
			node_normal = normalize(node_normal);
			node_to_normal_map[nidx] = node_normal;
			SLIC_INFO("normal "<< node_normal.m_x <<" " << node_normal.m_y << " " << node_normal.m_z);
		}

	}
*/

IndexType find_containing_zone(const axom::slam::IAMesh<3, 2>& mesh,
                               Point2D query_pt,
                               IndexType zone_i)
{
  SLIC_INFO("find_containing_zone " << zone_i);
  //assumes the point is definitely contained within a triangle
  //TODO fix that assumption ^
  SLIC_INFO("Query Pt " << query_pt);
  SLIC_ASSERT(mesh.dzn_rel[zone_i].size() == 3);  //should be a triangle...

  while(1)
  {
    SLIC_INFO("step into zone" << zone_i);
    IndexType n0 = mesh.dzn_rel[zone_i][0];
    IndexType n1 = mesh.dzn_rel[zone_i][1];
    IndexType n2 = mesh.dzn_rel[zone_i][2];

    Triangle2D tri(mesh.dcoord_map[n0], mesh.dcoord_map[n1], mesh.dcoord_map[n2]);

    axom::primal::Point3D bary_co = tri.barycentricCoords(query_pt);
    SLIC_INFO("bary_co " << bary_co);

    //Find the most negative
    IndexType i = 0;
    if(bary_co[0] > bary_co[1]) i = 1;
    if(bary_co[i] > bary_co[2]) i = 2;

    if(bary_co[i] >= 0)
    {                 //smaller than i -> outside of the triangle
      return zone_i;  //return if inside or on triangle
    }

    zone_i = mesh.dzz_rel[zone_i][2 - i];
  }
}

bool is_point_in_circle(Point2D p0, Point2D p1, Point2D p2, Point2D q)
{  //todo const
  double det = axom::numerics::determinant(1.0,
                                           p0[0],
                                           p0[1],
                                           p0[0] * p0[0] + p0[1] * p0[1],
                                           1.0,
                                           p1[0],
                                           p1[1],
                                           p1[0] * p1[0] + p1[1] * p1[1],
                                           1.0,
                                           p2[0],
                                           p2[1],
                                           p2[0] * p2[0] + p2[1] * p2[1],
                                           1.0,
                                           q[0],
                                           q[1],
                                           q[0] * q[0] + q[1] * q[1]);
  return det < 0;
}

/* Find the zones around a given zone whose delaunay circle contains the query point.
	 *
	 */
std::vector<IndexType> find_violating_zones(const axom::slam::IAMesh<3, 2>& mesh,
                                            Point2D query_pt,
                                            IndexType q_zone_i)
{
  //SLIC_INFO("find_violating_zones from zone " << q_zone_i );
  IndexType zone_i = q_zone_i;
  std::vector<IndexType> ret;
  std::list<IndexType> zone_list_to_check;
  zone_list_to_check.push_back(zone_i);

  std::set<IndexType> checked_zones;
  checked_zones.insert(zone_i);
  checked_zones.insert(-1);

  //starting from zone_i, which contains the point, try its neighbors, and the neighbor's neighbors
  //for point in circle test

  while(zone_list_to_check.size() > 0)
  {
    zone_i = zone_list_to_check.front();
    zone_list_to_check.pop_front();
    if(zone_i < 0) continue;
    IndexType n0 = mesh.dzn_rel[zone_i][0];
    IndexType n1 = mesh.dzn_rel[zone_i][1];
    IndexType n2 = mesh.dzn_rel[zone_i][2];
    Point2D p0 = mesh.dcoord_map[n0];
    Point2D p1 = mesh.dcoord_map[n1];
    Point2D p2 = mesh.dcoord_map[n2];
    bool is_in_circle = is_point_in_circle(p0, p1, p2, query_pt);
    if(is_in_circle)
    {
      bool is_repeat = false;
      if(!is_repeat)
      {
        ret.push_back(zone_i);
        for(int i = 0; i < (int)mesh.dzz_rel[zone_i].size(); i++)
        {
          IndexType nbr_zone = mesh.dzz_rel[zone_i][i];

          std::set<IndexType>::iterator it = checked_zones.find(nbr_zone);
          if(it == checked_zones.end())
          {
            zone_list_to_check.push_back(nbr_zone);
            checked_zones.insert(nbr_zone);
          }
        }
      }
    }
  }
  return ret;
}
}  // namespace slamTopologicalMesh

int printc = 1;
void print_out_mesh(axom::slam::IAMesh<3, 2> ia_mesh2)
{
  using namespace slamTopologicalMesh;
  SLIC_INFO("Printing out the resulting mesh:");
  std::cout << "i=" << printc++ << std::endl;
  std::cout << "nodes=[";
  for(int i = 0; i < ia_mesh2.dcoord_map.actual_size(); i++)
  {
    Point2D p = ia_mesh2.dcoord_map[i];
    std::cout << p[0] << "," << p[1] << ";" << std::endl;
  }
  std::cout << "];\n";
  std::cout << "tri=[";
  for(int i = 0; i < ia_mesh2.zone_set.actual_size(); i++)
  {
    if(ia_mesh2.zone_set.isValidEntry(i))
    {
      std::vector<int> nlist = ia_mesh2.getNodeInZone(i);
      std::cout << nlist[0] << "," << nlist[1] << "," << nlist[2] << ";"
                << std::endl;
    }
  }
  std::cout << "];\n";
}

int main(/*int argc, char** argv*/)
{
  using namespace slamTopologicalMesh;

  axom::slic::SimpleLogger logger;

  SLIC_INFO("Topological Mesh example");

  // Create an IA triangle mesh
  BasicTriMeshData basic_mesh_data;
  axom::slam::IAMesh<3, 3> ia_mesh(basic_mesh_data.points, basic_mesh_data.tri);

  //test zone->node function
  std::vector<IndexType> r;
  SLIC_INFO("list of nodes in zone 0");
  r = ia_mesh.getNodeInZone(0);
  for(unsigned int i = 0; i < r.size(); i++)
  {
    SLIC_INFO(r[i]);
  }

  SLIC_INFO("list of zones with node 0");
  r = ia_mesh.getZoneWithNode(0);
  for(unsigned int i = 0; i < r.size(); i++)
  {
    SLIC_INFO(r[i]);
  }

  //SLIC_INFO("Calculate normals for the nodes on the triangle mesh:");
  //calculate_normal(ia_mesh);

  //try it for tetrahedrons
  /*
	SLIC_INFO("Creating tetrahedron mesh");
	BasicTetMeshData basic_tet_mesh;
	axom::slam::IAMesh<4> ia_tetmesh(basic_tet_mesh.points, basic_tet_mesh.tet);

	SLIC_INFO("nodes in zone 0");
	r = ia_tetmesh.getNodeInZone(0);
	for(unsigned int i=0; i<r.size(); i++){
		SLIC_INFO(r[i]);
	}

	SLIC_INFO("zones with node 1");
	r = ia_tetmesh.getZoneWithNode(1);
	for(unsigned int i=0; i<r.size(); i++){
		SLIC_INFO(r[i]);
	}*/

  //SLIC_INFO("Find closest point triangles...");
  //find_node_closest_to_pt<3,2>(ia_mesh, axom::slam::Point(1,1,0));

  //find_node_closest_to_pt<4>(ia_tetmesh, axom::slam::Point(1,1,0));

  //randomly generate triangle mesh in 2d space
  Flat2DTriMeshData basic_tri_mesh_data;
  axom::slam::IAMesh<3, 2> ia_mesh2(basic_tri_mesh_data.points,
                                    basic_tri_mesh_data.tri);

  print_out_mesh(ia_mesh2);

  //fake adding a number of random points
  srand(1234);
  int num_points = 10;
  IndexType starting_zone = 0;
  for(int pt_i = 0; pt_i < num_points; pt_i++)
  {
    double x = (double)rand() / (double)RAND_MAX * 2.0 - 1.0;
    double y = (double)rand() / (double)RAND_MAX * 2.0 - 1.0;
    double pt_coord[2] = {x, y};
    Point2D new_pt(pt_coord, 2);

    //(assumes a mesh is already existing...)
    //find zone containing triangle
    IndexType zone_i = find_containing_zone(ia_mesh2, new_pt, starting_zone);

    //find the delaunay cavity...
    std::vector<IndexType> zones_to_remove =
      find_violating_zones(ia_mesh2, new_pt, zone_i);
    std::cout << "violating zones:";
    for(unsigned int i = 0; i < zones_to_remove.size(); i++)
      std::cout << zones_to_remove[i] << " ";
    std::cout << std::endl;

    //remove each zones from the mesh, keeping a list of cavity edges and their connected zones
    typedef std::pair<IndexType, IndexType> IndexPairType;
    typedef std::pair<IndexPairType, IndexType> CavityPairType;
    typedef std::map<IndexPairType, IndexType> CavityMapType;
    CavityMapType cavity_edges;
    for(unsigned int i = 0; i < zones_to_remove.size(); i++)
    {
      IndexType zone_i = zones_to_remove[i];
      for(int j = 0; j < (int)ia_mesh2.dzn_rel[zone_i].size(); j++)
      {
        IndexType z = zone_i;
        IndexType n1 = ia_mesh2.dzn_rel[z][(j * 2) % 3];
        IndexType n2 = ia_mesh2.dzn_rel[z][(1 + j * 2) % 3];
        CavityMapType::iterator iter = cavity_edges.find(IndexPairType(n2, n1));
        if(iter == cavity_edges.end())
        {
          cavity_edges.insert(CavityPairType(IndexPairType(n1, n2), zone_i));
        }
        else
        {  //already exists
          cavity_edges.erase(iter);
        }
      }
      ia_mesh2.removeZone(zone_i);
    }

    SLIC_INFO("After removing zones");
    print_out_mesh(ia_mesh2);
    std::cout << "pt = [" << pt_coord[0] << " " << pt_coord[1] << "];"
              << std::endl;

    //new triangles from the cavity edges

    IndexType new_pt_i = ia_mesh2.addNode(new_pt);
    for(CavityMapType::iterator iter = cavity_edges.begin();
        iter != cavity_edges.end();
        iter++)
    {
      IndexType n1 = iter->first.first;
      IndexType n2 = iter->first.second;
      //IndexType nbr_zone_i = iter->second;

      starting_zone = ia_mesh2.addZone(n1, n2, new_pt_i);
      SLIC_INFO("New cell: " << starting_zone << " with " << n1 << ", " << n2
                             << ", " << new_pt_i);
    }
    SLIC_INFO("After adding zones");
    print_out_mesh(ia_mesh2);
    std::cout << "pt = [" << pt_coord[0] << " " << pt_coord[1] << "];"
              << std::endl;
  }

  print_out_mesh(ia_mesh2);

  SLIC_INFO("Done!");

  return 0;
}
