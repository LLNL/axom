// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#include "runMIR.hpp"

int runMIR_seq(const conduit::Node &mesh,
               const conduit::Node &options,
               conduit::Node &result)
{
  return runMIR<axom::SEQ_EXEC>(mesh, options, result);
}
