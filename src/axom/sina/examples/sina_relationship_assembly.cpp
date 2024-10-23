// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Create IDs for both Task 22 and Run 1024
  axom::sina::ID task22 {"Task_22", axom::sina::IDType::Global};
  axom::sina::ID run1024 {"Run_1024", axom::sina::IDType::Global};

  // Create the relationship and print it out
  axom::sina::Relationship myRelationship {task22, "contains", run1024};
  std::cout << myRelationship.toNode().to_json() << std::endl;
}
