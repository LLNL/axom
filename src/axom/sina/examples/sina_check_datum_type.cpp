// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Define 3 different datums
  axom::sina::Datum myDatum {12.34};
  std::string value = "foobar";
  axom::sina::Datum myOtherDatum {value};
  std::vector<double> scalars = {1, 2, 20.0};
  axom::sina::Datum myArrayDatum {scalars};

  // Prints 0, corresponding to string
  std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(
                 myDatum.getType())
            << std::endl;

  // Prints 1, corresponding to scalar
  std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(
                 myOtherDatum.getType())
            << std::endl;

  // Prints 3, corresponding to scalar array
  std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(
                 myArrayDatum.getType())
            << std::endl;
}