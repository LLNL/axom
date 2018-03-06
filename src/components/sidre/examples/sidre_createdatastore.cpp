/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/* This example code contains snippets used in the Sidre Sphinx documentation.
 * They begin and end with comments
 *
 * first_example_create_start
 * first_example_create_end
 * first_example_chain_1
 * first_example_chain_2
 * first_example_access_start
 * first_example_access_end
 * serial_io_save_start
 * serial_io_save_end
 * tiny_create_start
 * tiny_create_end
 * blueprint_restructure_start
 * blueprint_restructure_end
 * blueprint_save_start
 * blueprint_save_end
 *
 * each prepended with an underscore.
 */

#include <cstring>

#include "sidre/sidre.hpp"

// "using" directives to simplify code
using namespace axom;
using namespace sidre;

DataStore * create_datastore(int * region) {

  // _first_example_create_start
  // Create Sidre datastore object and get root group
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

  // Create two attributes
  ds->createAttributeScalar("vis", 0);
  ds->createAttributeScalar("restart", 1);

  // Create group children of root group
  Group* state = root->createGroup("state");
  Group* nodes = root->createGroup("nodes");
  Group* fields = root->createGroup("fields");

  // Populate "state" group
  state->createViewScalar("cycle", 25);
  state->createViewScalar("time", 1.2562e-2);
  state->createViewString("name", "sample_20171206_a");

  int N = 16;
  int nodecount = N * N * N;
  int eltcount = (N-1) * (N-1) * (N-1);

  // Populate "nodes" group
  // 
  // "x", "y", and "z" are three views into a shared Sidre buffer object that
  // holds 3 * nodecount doubles.  These views might describe the location of
  // each node in a 16 x 16 x 16 hexahedron mesh.  Each view is described by
  // number of elements, offset, and stride into that data.
  // _first_example_chain_1
  Buffer* buff = ds->createBuffer(sidre::DOUBLE_ID, 3*nodecount)->allocate();
  nodes->createView("x", buff)->apply(sidre::DOUBLE_ID, nodecount, 0, 3);
  // _first_example_chain_2
  nodes->createView("y", buff)->apply(sidre::DOUBLE_ID, nodecount, 1, 3);
  nodes->createView("z", buff)->apply(sidre::DOUBLE_ID, nodecount, 2, 3);

  // Populate "fields" group
  //
  // "temp" is a view into a buffer that is not shared with another View.
  // In this case, the data Buffer is allocated directly through the View
  // object.  Likewise with "rho."  Both Views have the default offset (0)
  // and stride (1).  These Views could point to data associated with
  // each of the 15 x 15 x 15 hexahedron elements defined by the nodes above.
  View * temp =
    fields->createViewAndAllocate("temp", sidre::DOUBLE_ID, eltcount);
  View * rho =
    fields->createViewAndAllocate("rho", sidre::DOUBLE_ID, eltcount);

  // Explicitly set values for the "vis" Attribute on the "temp" and "rho"
  // buffers.
  temp->setAttributeScalar("vis", 1);
  rho->setAttributeScalar("vis", 1);

  // The "fields" Group also contains a child Group "ext" which holds a pointer
  // to an externally owned integer array.  Although Sidre does not own the
  // data, the data can still be described to Sidre.
  Group* ext = fields->createGroup("ext");
  // int * region has been passed in as a function argument.  As with "temp"
  // and "rho", view "region" has default offset and stride.
  ext->createView("region", region)->apply(sidre::INT_ID, eltcount);
  // _first_example_create_end

  return ds;
}

void access_datastore(DataStore * ds) {
  // _first_example_access_start
  // Retrieve Group pointers
  Group * root = ds->getRoot();
  Group * state = root->getGroup("state");
  Group * nodes = root->getGroup("nodes");
  Group * fields = root->getGroup("fields");

  // Access items in "state" group
  int cycle = state->getView("cycle")->getScalar();
  double time = state->getView("time")->getScalar();
  const char* name = state->getView("name")->getString();

  // Access some items in "nodes" and "fields" groups
  double* y = nodes->getView("y")->getArray();
  int ystride = nodes->getView("y")->getStride();
  double* temp = fields->getView("temp")->getArray();
  int* region = fields->getView("ext/region")->getArray();

  // Nudge the 3rd node, adjust temp and region of the 3rd element
  y[2 * ystride] += 0.0032;
  temp[2] *= 1.0021;
  region[2] = 6;
  // _first_example_access_end

  // Deal with unused variables
  AXOM_DEBUG_VAR(cycle);
  AXOM_DEBUG_VAR(time);
  AXOM_DEBUG_VAR(name);
  AXOM_DEBUG_VAR(y);
  AXOM_DEBUG_VAR(temp);
  AXOM_DEBUG_VAR(region);
}

DataStore * create_tiny_datastore() {
  // _tiny_create_start
  DataStore *ds = new DataStore();

  int nodecount = 12;
  int elementcount = 2;

  // Create views and buffers to hold node positions and field values
  Group * nodes = ds->getRoot()->createGroup("nodes");
  View * xs = nodes->createViewAndAllocate("xs", sidre::DOUBLE_ID, nodecount);
  View * ys = nodes->createViewAndAllocate("ys", sidre::DOUBLE_ID, nodecount);
  View * zs = nodes->createViewAndAllocate("zs", sidre::DOUBLE_ID, nodecount);

  Group * fields = ds->getRoot()->createGroup("fields");
  View * nodefield =
    fields->createViewAndAllocate("nodefield", sidre::INT_ID, nodecount);
  View * eltfield =
    fields->createViewAndAllocate("eltfield", sidre::DOUBLE_ID, elementcount);

  // Set node position for two adjacent hexahedrons
  double * xptr = xs->getArray();
  double * yptr = ys->getArray();
  double * zptr = zs->getArray();
  for (int pos = 0; pos < nodecount; ++pos) {
    xptr[pos] = (pos / 2) % 2;
    yptr[pos] = ((pos + 1) / 2) % 2;
    zptr[pos] = zpos / 4;
  }

  // Assign a value to the node field
  int * nf = nodefield->getArray();
  for (int pos = 0; pos < nodecount; ++pos) {
    nf[pos] = 8 - pos;
  }
  // and to the element field.
  double * ef = eltfield->getArray();
  // There are only two.
  ef[0] = 2.65;
  ef[1] = 1.96;

  return ds;
  // _tiny_create_end
}

void save_as_blueprint(DataStore * ds) {
  // _blueprint_restructure_start
  // Conduit needs a specific hierarchy.
  // We'll make a new DataStore with that hierarchy, pointing at the
  // application's data.
  DataStore *cds = new DataStore();
  
  // The Conduit specifies top-level groups:
  Group * coords = cds->getRoot()->createGroup("coordsets/coords");
  Group * topos = cds->getRoot()->createGroup("topologies");
  // no material sets in this example
  Group * fields = cds->getRoot()->createGroup("fields");
  Group * adj = cds->getRoot()->createGroup("adjsets");

  // Set up the coordinates as Mesh Blueprint requires
  coords->createViewString("type", "explicit");
  // We use prior knowledge of the layout of the original datastore
  View * origv = ds->getRoot()->getView("nodes/xs");
  View * conduitval = coords->createGroup("values");
  conduitval->createView("x", origv->getTypeID(),
                         origv->getNumElements(),
                         origv->getArray());
  origv = ds->getRoot()->getView("nodes/ys");
  conduitval->createView("y", origv->getTypeID(),
                         origv->getNumElements(),
                         origv->getArray());
  origv = ds->getRoot()->getView("nodes/zs");
  conduitval->createView("z", origv->getTypeID(),
                         origv->getNumElements(),
                         origv->getArray());
  // _blueprint_restructure_end
}

void serial_save_datastore_and_load_copy_lower(DataStore *ds) {
  // _serial_io_save_start
  // Save the data store to a file, using the default sidre_hdf5 protocol,
  // saving all Views
  ds->getRoot()->save("example.hdf5");
  // Delete the data hierarchy under the root, then load it from the file
  ds->getRoot()->load("example.hdf5");
  Group * additional = ds->getRoot()->createGroup("additional");
  additional->createGroup("yetanother");
  // Load another copy of the data store into the "additional" group
  // without first clearing all its contents
  additional->load("example.hdf5", "sidre_hdf5", true);
  // _serial_io_save_end
}

int main(int argc, char ** argv) {

  // Deal with unused variables
  AXOM_DEBUG_VAR(argc);
  AXOM_DEBUG_VAR(argv);

  int region[3375];

  DataStore * ds = create_datastore(region);
  access_datastore(ds);

  return 0;
}
