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
#include "axom/slic.hpp"

#include "axom/fmt.hpp"
#include "axom/CLI11.hpp"

#include <cmath>
#include <map>
#include <string>

using HexahedronType = typename axom::primal::Hexahedron<double>;
using PointType = typename axom::primal::Point<double, 3>;
using TetrahedronType = typename axom::primal::Tetrahedron<double, 3>;

/// Choose runtime policy for RAJA
enum RuntimePolicy
{
  seq = 0,
  omp = 1,
  #ifdef AXOM_USE_CUDA
  cuda = 2,
  #endif
  #ifdef AXOM_USE_HIP
  hip = 3
  #endif
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
      {"seq", seq}
    #if defined(AXOM_USE_RAJA)
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

// Helper function to generate a hexahedron cube with given corner point and
// edge length
HexahedronType generateCube(const PointType& point, double length)
{
  return HexahedronType(
    point,
    PointType::make_point(point[0], point[1], point[2] + length),
    PointType::make_point(point[0] + length, point[1], point[2] + length),
    PointType::make_point(point[0] + length, point[1], point[2]),
    PointType::make_point(point[0], point[1] + length, point[2]),
    PointType::make_point(point[0], point[1] + length, point[2] + length),
    PointType::make_point(point[0] + length, point[1] + length, point[2] + length),
    PointType::make_point(point[0] + length, point[1] + length, point[2]));
}

// Function to check intersection volumes of generated hexahedra and tetrahedra
template <typename ExecSpace>
void check_intersection_volumes(const Input& params)
{
  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(
      "Running intersection volume check in execution Space: {}",
      axom::execution_space<ExecSpace>::name())));

  // Save current/default allocator
  const int current_allocator = axom::getDefaultAllocatorID();

  // Set new default allocator
  axom::setDefaultAllocator(axom::execution_space<ExecSpace>::allocatorID());

  // Generate hexahedra subdividing the unit cube with corner points (-1,-1,-1) and (1,1,1)
  // and corresponding bounding boxes
  int const HEX_LEVEL = 2;
  int hex_index = 0;
  int const NUM_HEXES = HEX_LEVEL * HEX_LEVEL * HEX_LEVEL;
  HexahedronType* hexes = axom::allocate<HexahedronType>(NUM_HEXES);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Generating {} hexahedra with hexahedra level set to {}",
                      NUM_HEXES,
                      HEX_LEVEL)));

  for(int i = 0; i < HEX_LEVEL; i++)
  {
    for(int j = 0; j < HEX_LEVEL; j++)
    {
      for(int k = 0; k < HEX_LEVEL; k++)
      {
        double edge_length = 2.0 / HEX_LEVEL;
        hexes[hex_index] =
          generateCube(PointType::make_point(edge_length * i - 1,
                                             edge_length * j - 1,
                                             edge_length * k - 1),
                       edge_length);
        SLIC_INFO("Hex " << hex_index << " volume is "
                         << hexes[hex_index].volume() << "\n");
        SLIC_INFO("Hex " << hex_index << " is \n" << hexes[hex_index]);
        hex_index++;
      }
    }
  }

  // Generate tetrahedra from unit sphere with center (0,0,0) and corresponding bounding boxes
  int const TET_LEVEL = 1;
  int tet_index = 0;
  int const NUM_TETS = 4 * std::pow(2, TET_LEVEL);
  TetrahedronType* tets = axom::allocate<TetrahedronType>(NUM_TETS);
  double step = 1.0 / (TET_LEVEL * 2);

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(
      "Generating {} tetrahedra with tetrahedra level set to {}",
      NUM_TETS,
      TET_LEVEL)));

  for(int i = 0; i < TET_LEVEL * 2; i++)
  {
    for(int j = 0; j <= 1; j++)
    {
      for(int k = 0; k <= 1; k++)
      {
        // Tetrahedron coordinates consist of two points on the sphere where y = 0,
        // the origin, and one of the "poles" of the unit sphere, either (0,1,0) or
        // (0, -1, 0)

        double pole_sign = j ? 1 : -1;
        double z_sign = k ? 1 : -1;

        double x1 = axom::utilities::lerp<double>(-1.0, 1.0, step * i);
        double z1 = std::sqrt(1 - (x1 * x1)) * z_sign;
        double x2 = axom::utilities::lerp<double>(-1.0, 1.0, step * (i + 1));
        double z2 = std::sqrt(1 - (x2 * x2)) * z_sign;

        tets[tet_index] = TetrahedronType(PointType::make_point(x1, 0, z1),
                                          PointType::zero(),
                                          PointType::make_point(x2, 0, z2),
                                          PointType::make_point(0, pole_sign, 0));
        SLIC_INFO("Tet " << tet_index << " volume is "
                         << tets[tet_index].volume() << "\n");
        SLIC_INFO("Tet " << tet_index << " is \n" << tets[tet_index]);
        tet_index++;
      }
    }
  }

  // Calculate expected sum of all tetrahedra volume
  double total_tet_vol = 0.0;
  for(int i = 0; i < NUM_TETS; i++)
  {
    total_tet_vol += tets[i].volume();
  }

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format("Total volume of all tetrahedra is {}", total_tet_vol)));

  // May not need bounding boxes, if our algorithm is actually good,
  // and we can just do a clip() for each pair...

  double total_intersect_vol = 0.0;

  // Calculate intersection volume for each hexahedra and tetrahedra pair.
  for(int i = 0; i < NUM_HEXES; i++)
  {
    for(int j = 0; j < NUM_TETS; j++)
    {
      total_intersect_vol += intersection_volume(hexes[i], tets[j]);
      SLIC_INFO("Volume " << i << " " << j << " is "
                          << intersection_volume(hexes[i], tets[j]) << "\n");
    }
  }

  SLIC_INFO(axom::fmt::format(
    "{:-^80}",
    axom::fmt::format(
      "Total intersect volume between all hexahedra and tetrahedra is {}",
      total_intersect_vol)));

  // Verify sum of intersection volumes is same as sum of all tetradra volumes

  // Reset default allocator
  axom::setDefaultAllocator(current_allocator);
}

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

  switch(params.policy)
  {
  case RuntimePolicy::seq:
    check_intersection_volumes<axom::SEQ_EXEC>(params);
    break;
#if defined(AXOM_USE_RAJA)
  #ifdef AXOM_USE_OPENMP
  case RuntimePolicy::omp:
    check_intersection_volumes<axom::OMP_EXEC>(params);
    break;
  #endif
  #ifdef AXOM_USE_CUDA
  case RuntimePolicy::cuda:
    check_intersection_volumes<axom::CUDA_EXEC<256>>(params);
    break;
  #endif
  #ifdef AXOM_USE_HIP
  case RuntimePolicy::hip:
    check_intersection_volumes<axom::HIP_EXEC<256>>(params);
    break;
  #endif
#endif
  }
  return 0;
}