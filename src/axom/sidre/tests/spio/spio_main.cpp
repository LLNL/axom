// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "gtest/gtest.h"

#include "mpi.h"

#include "axom/config.hpp"  // for compile-time definitions

namespace
{
#ifdef AXOM_USE_HDF5
const std::string PROTOCOL = "sidre_hdf5";
#else
const std::string PROTOCOL = "sidre_json";
#endif
const std::string ROOT_EXT = ".root";
}  // namespace

#include "spio_basic.hpp"
#include "spio_parallel.hpp"
#include "spio_serial.hpp"

#include "axom/slic/core/SimpleLogger.hpp"
using axom::slic::SimpleLogger;

int main(int argc, char* argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  SimpleLogger logger;  // create & initialize test logger,

  MPI_Init(&argc, &argv);
  result = RUN_ALL_TESTS();
  MPI_Finalize();

  return result;
}
