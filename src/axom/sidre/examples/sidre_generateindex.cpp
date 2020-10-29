// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*! \file sidre_createdatastore.cpp
 *  \brief This example code is a basic demonstration of how to use the
 *  Sidre classes to construct a hierarchical data store and save it to
 *  a file.  (This example does not use the parallel I/O facility.)
 *  This also shows how to use Sidre with Conduit to generate and save a
 *  data store using the Mesh Blueprint, for easy data exchange and
 *  visualization.
 */

/*
 * Summary of the functions called by main this example:
 *
 * create_datastore           -- Builds a DataStore containing various Groups,
 *                               Views, Buffers, and Attributes
 *
 * create_tiny_datastore      -- Builds a DataStore that holds a small 3D mesh
 *                               conforming to the conduit mesh blueprint
 *
 * access_datastore           -- Exercises access methods that get data from
 *                               the DataStore created by create_datastore
 *
 * save_as_blueprint          -- Generates a blueprint index and writes data
 *                               to file using direct calls to conduit
 *
 * generate_blueprint         -- Generates a blueprint index and writes data
 *                               to file using Sidre DataStore and Group calls
 *
 * generate_blueprint_to_path -- Similar to above, but uses path arguments to
 *                               find data deeper in Sidre hierarchy
 *
 * generate_spio_blueprint    -- Uses SPIO/IOManager calls to write data to
 *                               file and generate blueprint index
 *
 * generate_blueprint_to_path -- Similar to above, but uses path arguments to
 *                               find data deeper in Sidre hierarchy
 */

/* This example code contains snippets used in the Sidre Sphinx documentation.
 * They begin and end with comments
 *
 * first_example_creategroups_start
 * first_example_creategroups_end
 * first_example_state_start
 * first_example_state_end
 * first_example_nodes_start
 * first_example_nodes_end
 * first_example_fields_start
 * first_example_fields_end
 * first_example_access_start
 * first_example_access_end
 * serial_io_save_start
 * serial_io_save_end
 * tiny_create_start
 * tiny_create_end
 * blueprint_restructure_save_start
 * blueprint_restructure_save_end
 * blueprint_save_start
 * blueprint_save_end
 *
 * each prepended with an underscore.
 */

// Axom headers
#include "axom/sidre.hpp"

// Conduit headers
#include "conduit.hpp"
#include "conduit_blueprint.hpp"
#include "conduit_relay.hpp"

#ifdef AXOM_USE_MPI
  #include "conduit_blueprint_mpi.hpp"
#endif

//fmt format header
#include "fmt/fmt.hpp"

// C++ headers
#include <cstring>

// "using" directives to simplify code
using namespace axom;
using namespace sidre;

DataStore* create_tiny_datastore()
{
  // _tiny_create_start
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
  // _tiny_create_end
}

void setup_blueprint_coords(DataStore* ds, Group* coords)
{
  // _blueprint_restructure_coords_start
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
  // _blueprint_restructure_coords_end
}

void setup_blueprint_external_coords(DataStore* ds, Group* coords)
{
  // _blueprint_external_coords_start
  // Set up the coordinates as Mesh Blueprint requires
  coords->createViewString("type", "explicit");
  // We use prior knowledge of the layout of the original datastore
  View* origv = ds->getRoot()->getView("nodes/xs");
  Group* conduitval = coords->createGroup("values");
  conduitval->createView("x",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  origv = ds->getRoot()->getView("nodes/ys");
  conduitval->createView("y",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  origv = ds->getRoot()->getView("nodes/zs");
  conduitval->createView("z",
                         sidre::DOUBLE_ID,
                         origv->getNumElements(),
                         static_cast<double*>(origv->getArray()));
  // _blueprint_external_coords_end
}

void setup_blueprint_topos(DataStore* ds, Group* topos)
{
  // _blueprint_restructure_topo_start
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
  // _blueprint_restructure_topo_end

  // Deal with unused variable
  AXOM_DEBUG_VAR(ds);
}

void setup_blueprint_fields(DataStore* ds, Group* fields)
{
  // _blueprint_restructure_field_start
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
  // _blueprint_restructure_field_end
}

void setup_blueprint_external_fields(DataStore* ds, Group* fields)
{
  // _blueprint_external_field_start
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
                        static_cast<int*>(origv->getArray()));

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
                       static_cast<double*>(origv->getArray()));
  // _blueprint_external_field_end
}

#ifdef AXOM_USE_MPI
void generate_spio_blueprint(DataStore* ds, bool dense)
{

  int my_rank;
  int comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

//  std::string domain_name = fmt:sprintf("%s_%06d", "domain", my_rank);
  std::string domain_name = "domain";
  std::string holder_name = "domain_data";
  std::string domain_location = holder_name + "/" + domain_name;
  std::string mesh_name = "mesh";

  Group* holder = ds->getRoot()->createGroup(holder_name);

  // dense places domains on all ranks when true. When not true, it
  // leaves off rank zero to show that the index can still be generated
  // with some ranks empty.
  if (dense || my_rank > 0 || comm_size == 0)
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

  // _blueprint_generate_spio_start
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
    if (dense)
    {
      output_name = "bpdense";
    }
    else
    {
      output_name = "bpsparse";
    }
    std::string bp_rootfile = output_name + ".root";

    writer.write(ds->getRoot()->getGroup("domain_data"), 1, output_name, protocol);

    writer.writeBlueprintIndexToRootFile(ds, domain_location, bp_rootfile, mesh_name);
  }
  // _blueprint_generate_spio_end
}

void generate_spio_blueprint_to_path(DataStore* ds)
{
  // _blueprint_spio_path_start
  std::string domain_name = "domain";
  std::string domain_location = "domain_data/level/domains/" + domain_name;
  std::string mesh_name = "spiopathmesh";
  std::string domain_mesh = domain_location + "/" + mesh_name;

  Group* mroot = ds->getRoot()->createGroup(domain_location);
  Group* coords = mroot->createGroup(mesh_name + "/coordsets/coords");
  Group* topos = mroot->createGroup(mesh_name + "/topologies");
  // no material sets in this example
  Group* fields = mroot->createGroup(mesh_name + "/fields");
  // no adjacency sets in this (single-domain) example
  // _blueprint_spio_path_end

  setup_blueprint_external_coords(ds, coords);

  setup_blueprint_topos(ds, topos);

  setup_blueprint_external_fields(ds, fields);

  // _blueprint_generate_spio_path_start
  IOManager writer(MPI_COMM_WORLD);

  conduit::Node info, mesh_node, root_node;
  ds->getRoot()->createNativeLayout(mesh_node);
  std::string bp_protocol = "mesh";
  if(conduit::blueprint::mpi::verify(bp_protocol,
                                     mesh_node[domain_mesh],
                                     info,
                                     MPI_COMM_WORLD))
  {
    std::string bp_rootfile("pathbpspio.root");
  #if defined(AXOM_USE_HDF5)
    std::string protocol = "sidre_hdf5";
  #else
    std::string protocol = "sidre_json";
  #endif

    writer.write(ds->getRoot()->getGroup("domain_data"), 1, "pathbpspio", protocol);

    writer.writeBlueprintIndexToRootFile(ds,
                                         domain_mesh,
                                         bp_rootfile,
                                         "level/domains/domain/" + mesh_name);
  }
  // _blueprint_generate_spio_path_end
}
#endif

int main(int argc, char** argv)
{
  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);
#ifdef AXOM_USE_MPI
  int num_ranks = 1;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  DataStore* sds = create_tiny_datastore();
  DataStore* spds = create_tiny_datastore();
  generate_spio_blueprint(sds, true);   //dense
  generate_spio_blueprint(spds, false); //sparse
//  generate_spio_blueprint_to_path(spds);
  MPI_Finalize();
#endif

  return 0;
}
