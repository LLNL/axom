// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file sidre_generateindex.cpp
 *  \brief This example code is a basic demonstration of the generation
 *  of a Sidre DataStore containing a mesh conforming to the Mesh Blueprint,
 *  including the automatic generation of a Blueprint index and saving it
 *  to a file.
 */

/*
 * Summary of the functions called by main this example:
 *
 * create_tiny_datastore      -- Builds a DataStore that holds a small 3D mesh
 *                               conforming to the conduit mesh blueprint
 *
 * generate_spio_blueprint    -- Uses SPIO/IOManager calls to write data to
 *                               file and generate blueprint index
 */

// Axom headers
#include "axom/core.hpp"
#include "axom/sidre.hpp"

// Conduit headers
#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
#endif

// fmt format header
#include "fmt/fmt.hpp"

// C++ headers
#include <cstring>

// "using" directives to simplify code
using namespace axom;
using namespace sidre;

DataStore* create_tiny_datastore()
{
  DataStore* ds = new DataStore();

  int nodecount = 12;
  int elementcount = 2;

  // Create views and buffers to hold node positions and field values
  Group* nodes = ds->getRoot()->createGroup("nodes");
  View* xs = nodes->createViewAndAllocate("xs", sidre::DOUBLE_ID, nodecount);
  View* ys = nodes->createViewAndAllocate("ys", sidre::DOUBLE_ID, nodecount);
  View* zs = nodes->createViewAndAllocate("zs", sidre::DOUBLE_ID, nodecount);

  Group* fields = ds->getRoot()->createGroup("fields");
  View* nodefield =
    fields->createViewAndAllocate("nodefield", sidre::INT_ID, nodecount);
  View* eltfield =
    fields->createViewAndAllocate("eltfield", sidre::DOUBLE_ID, elementcount);

  // Set node position for two adjacent hexahedrons
  double* xptr = xs->getArray();
  double* yptr = ys->getArray();
  double* zptr = zs->getArray();
  for(int pos = 0; pos < nodecount; ++pos)
  {
    xptr[pos] = ((pos + 1) / 2) % 2;
    yptr[pos] = (pos / 2) % 2;
    zptr[pos] = pos / 4;
  }

  // Assign a value to the node field
  int* nf = nodefield->getArray();
  for(int pos = 0; pos < nodecount; ++pos)
  {
    nf[pos] = static_cast<int>(xptr[pos] + yptr[pos] + zptr[pos]);
  }
  // and to the element field.
  double* ef = eltfield->getArray();
  // There are only two elements.
  ef[0] = 2.65;
  ef[1] = 1.96;

  return ds;
}

void setup_blueprint_coords(DataStore* ds, Group* coords)
{
  // Set up the coordinates as Mesh Blueprint requires
  coords->createViewString("type", "explicit");
  // We use prior knowledge of the layout of the original datastore
  View* origv = ds->getRoot()->getView("nodes/xs");
  Group* conduitval = coords->createGroup("values");
  conduitval->createView("x",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         origv->getBuffer());
  origv = ds->getRoot()->getView("nodes/ys");
  conduitval->createView("y",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         origv->getBuffer());
  origv = ds->getRoot()->getView("nodes/zs");
  conduitval->createView("z",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         origv->getBuffer());
}

void setup_blueprint_topos(DataStore* ds, Group* topos)
{
  // Sew the nodes together into the two hexahedra, using prior knowledge.
  Group* connmesh = topos->createGroup("mesh");
  connmesh->createViewString("type", "unstructured");
  connmesh->createViewString("coordset", "coords");
  Group* elts = connmesh->createGroup("elements");
  elts->createViewString("shape", "hex");

  // We have two eight-node hex elements, so we need 2 * 8 = 16 ints.
  View* connectivity =
    elts->createViewAndAllocate("connectivity", sidre::INT_ID, 16);

  // The Mesh Blueprint connectivity array for a hexahedron lists four nodes on
  // one face arranged by right-hand rule to indicate a normal pointing into
  // the element, then the four nodes of the opposite face arranged to point
  // the normal the same way (out of the element).  This is the same as for
  // a VTK_HEXAHEDRON.  See
  // https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf.

  int* c = connectivity->getArray();

  // First hex.  In this example, the Blueprint node ordering matches the
  // dataset layout.  This is fortuitous but not required.
  c[0] = 0;
  c[1] = 1;
  c[2] = 2;
  c[3] = 3;
  c[4] = 4;
  c[5] = 5;
  c[6] = 6;
  c[7] = 7;

  // Second and last hex
  c[8] = 4;
  c[9] = 5;
  c[10] = 6;
  c[11] = 7;
  c[12] = 8;
  c[13] = 9;
  c[14] = 10;
  c[15] = 11;

  // Deal with unused variable
  AXOM_UNUSED_VAR(ds);
}

void setup_blueprint_fields(DataStore* ds, Group* fields)
{
  // Set up the node-centered field
  // Get the original data
  View* origv = ds->getRoot()->getView("fields/nodefield");
  Group* nodefield = fields->createGroup("nodefield");
  nodefield->createViewString("association", "vertex");
  nodefield->createViewString("type", "scalar");
  nodefield->createViewString("topology", "mesh");
  nodefield->createView("values",
                        sidre::INT_ID,
                        origv->getNumElements(),
                        origv->getBuffer());

  // Set up the element-centered field
  // Get the original data
  origv = ds->getRoot()->getView("fields/eltfield");
  Group* eltfield = fields->createGroup("eltfield");
  eltfield->createViewString("association", "element");
  eltfield->createViewString("type", "scalar");
  eltfield->createViewString("topology", "mesh");
  eltfield->createView("values",
                       sidre::DOUBLE_ID,
                       origv->getNumElements(),
                       origv->getBuffer());
}

#ifdef AXOM_USE_MPI
void generate_spio_blueprint(DataStore* ds, bool dense)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  std::string domain_name = "domain";
  std::string holder_name = "domain_data";
  std::string domain_location = holder_name + "/" + domain_name;
  std::string mesh_name = "mesh";

  Group* holder = ds->getRoot()->createGroup(holder_name);

  // dense places domains on all ranks when true. When not true, it
  // leaves off rank zero to show that the index can still be generated
  // with some ranks empty.
  if(dense || my_rank > 0 || comm_size == 1)
  {
    Group* mroot = holder->createGroup(domain_name);
    Group* coords = mroot->createGroup("coordsets/coords");
    Group* topos = mroot->createGroup("topologies");
    // no material sets in this example
    Group* fields = mroot->createGroup("fields");
    // no adjacency sets in this (single-domain) example

    setup_blueprint_coords(ds, coords);

    setup_blueprint_topos(ds, topos);

    setup_blueprint_fields(ds, fields);
  }

  IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol,
                                     mesh_node[domain_location],
                                     info,
                                     MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    std::string output_name;
    if(dense)
    {
      output_name = "bpdense";
    }
    else
    {
      output_name = "bpsparse";
    }

    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"), 1, output_name, protocol);

    writer.writeBlueprintIndexToRootFile(ds,
                                         domain_location,
                                         bp_rootfile,
                                         mesh_name);
  }
}

void generate_multidomain_blueprint(DataStore* ds,
                                    const std::string& filename,
                                    int num_files)
{
  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  // 3 domains on even ranks, 2 domains on odd ranks
  int64_t domain_begin = (5 * (my_rank / 2)) + (3 * (my_rank % 2));
  int64_t domain_end = domain_begin + (3 - (my_rank % 2));

  std::string holder_name = "domain_data";
  std::string mesh_name = "mesh";
  Group* holder = ds->getRoot()->createGroup(holder_name);

  std::string domain_pattern = "domain_{domain:06d}";
  ds->getRoot()->createViewString("domain_pattern", domain_pattern);

  for(int64_t i = domain_begin; i < domain_end; ++i)
  {
    std::string domain_name = fmt::format("domain_{:06d}", i);

    Group* mroot = holder->createGroup(domain_name);
    Group* coords = mroot->createGroup("coordsets/coords");
    Group* topos = mroot->createGroup("topologies");
    // no material sets in this example
    Group* fields = mroot->createGroup("fields");
    // no adjacency sets in this example
    mroot->createViewScalar("state/domain_id", i);

    setup_blueprint_coords(ds, coords);

    setup_blueprint_topos(ds, topos);

    setup_blueprint_fields(ds, fields);
  }

  IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol,
                                     mesh_node[holder_name],
                                     info,
                                     MPI_COMM_WORLD))
  {
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    std::string output_name = filename;

    if(comm_size > 1)
    {
      output_name = output_name + "_par";
    }

    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"),
                 std::min(num_files, comm_size),
                 output_name,
                 protocol);

    writer.writeBlueprintIndexToRootFile(ds, holder_name, bp_rootfile, mesh_name);
  }
}
#endif

int main(int argc, char** argv)
{
  AXOM_UNUSED_VAR(argc);
  AXOM_UNUSED_VAR(argv);
#ifdef AXOM_USE_MPI
  int num_ranks = 1;

  slic::SimpleLogger logger;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  DataStore* sds = create_tiny_datastore();
  DataStore* spds = create_tiny_datastore();
  generate_spio_blueprint(sds, true);    //dense
  generate_spio_blueprint(spds, false);  //sparse
  spds->getRoot()->destroyGroups();
  spds->getRoot()->destroyViews();
  spds = create_tiny_datastore();
  generate_multidomain_blueprint(spds, "multi1", 1);
  if(num_ranks > 1)
  {
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();
    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds, "multi2", 2);
    spds->getRoot()->destroyGroups();
    spds->getRoot()->destroyViews();
    spds = create_tiny_datastore();
    generate_multidomain_blueprint(spds, "multi_all", num_ranks);
  }
  MPI_Finalize();
#endif

  return 0;
}
