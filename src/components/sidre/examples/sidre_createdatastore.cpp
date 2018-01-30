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

#include <cstring>

#include "sidre/sidre.hpp"

// "using" directives to simplify code
using namespace axom;
using namespace sidre;

DataStore * create_datastore(int * region) {
  // first_example_create_start

  // Create Sidre datastore object and get root group
  DataStore* ds = new DataStore();
  Group* root = ds->getRoot();

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
  // first_example_chain_1
  Buffer* buff = ds->createBuffer(sidre::DOUBLE_ID, 3*nodecount)->allocate();
  // first_example_chain_2
  nodes->createView("x", buff)->apply(sidre::DOUBLE_ID, nodecount, 0, 3);
  // first_example_chain_3
  nodes->createView("y", buff)->apply(sidre::DOUBLE_ID, nodecount, 1, 3);
  nodes->createView("z", buff)->apply(sidre::DOUBLE_ID, nodecount, 2, 3);

  // Populate "fields" group
  //
  // "temp" is a view into a buffer that is not shared with another View.
  // In this case, the data Buffer can be allocated directly through the View
  // object.  Likewise with "rho."  These might contain data associated with
  // each of the 15 x 15 x 15 hexahedron elements defined by the nodes above.
  fields->createViewAndAllocate("temp", sidre::DOUBLE_ID, eltcount);
  fields->createViewAndAllocate("rho", sidre::DOUBLE_ID, eltcount);

  // The "fields" Group also contains a child Group "ext" which holds a pointer
  // to an externally owned integer array.  Although Sidre does not own the
  // data, the data can still be described to Sidre.
  Group* ext = fields->createGroup("ext");
  // int * region has been passed in as a function argument.
  ext->createView("region", region)->apply(sidre::INT_ID, eltcount);
  // first_example_create_end

  return ds;
}

void access_datastore(DataStore * ds) {
  // first_example_access_start
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
  double* temp = fields->getView("temp")->getArray();
  int* region = fields->getView("ext/region")->getArray();
  // first_example_access_end
}

int main(int argc, char ** argv) {
  int region[3375];

  DataStore * ds = create_datastore(region);
  access_datastore(ds);

  return 0;
}
