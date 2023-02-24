// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 * \file hex_tet_volume.cpp
 *
 * Example that demonstrates use of Primal's intersection_volume operator to find
 * the volume of intersection between hexahedra and tetrahedra.  Supports host
 * and device execution using RAJA.
 */

// Idea: have two input values to control how many hexes to make to represent the hex mesh
// and how many tets to make to revolve around the sphere
// and one input value for where to run this operation host or device
// Goal is to check if clip volumes sum is same as total tet volumes, within a tolerance.
// Not using spin or quest to have this located in primal.

#include "axom/primal.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

/// Choose runtime policy for RAJA
enum RuntimePolicy
{
  seq = 0,
  omp = 1,
  cuda = 2,
  hip = 3
};

struct Input
{
public:
  RuntimePolicy policy {seq};

  void parse(int argc, char** argv, axom::CLI::App& app)
  {
    std::string pol_str =
      "Sets execution space of intersection_volume operator.";
#if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
    pol_str += "\nSet to \'seq\' or 0 to use the RAJA sequential policy.";
  #ifdef AXOM_USE_OPENMP
    pol_str += "\nSet to \'omp\' or 1 to use the RAJA OpenMP policy.";
  #endif
  #ifdef AXOM_USE_CUDA
    pol_str += "\nSet to \'cuda\' or 2 to use the RAJA CUDA policy.";
  #endif
  #ifdef AXOM_USE_HIP
    pol_str += "\nSet to \'hip\' or 3 to use the RAJA HIP policy.";
  #endif
#endif

    app.add_option("-p, --policy", this->policy, pol_str)
      ->capture_default_str()
      ->transform(axom::CLI::CheckedTransformer(validExecPolicies));

    app.get_formatter()->column_width(50);

    app.parse(argc, argv);
  }

private:
  // clang-format off
  const std::map<std::string, RuntimePolicy> validExecPolicies{
    #if defined(AXOM_USE_RAJA) && defined(AXOM_USE_UMPIRE)
      {"seq", seq}
      #ifdef AXOM_USE_OPENMP
    , {"omp", omp}
      #endif
      #ifdef AXOM_USE_CUDA
    , {"cuda", cuda}
      #endif
      #ifdef AXOM_USE_HIP
    , {"hip", hip}
      #endif
    #endif
  };
  // clang-format on
};

int main(int argc, char** argv)
{
  // Initialize the SLIC logger
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  // Set up and parse command line arguments
  Input params;
  axom::CLI::App app {
    "Example of intersection volume between hexahedra and tetrahedra"};

  try
  {
    params.parse(argc, argv, app);
  }
  catch(const axom::CLI::ParseError& e)
  {
    return app.exit(e);
  }

  return 0;
}