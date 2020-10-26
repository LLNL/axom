// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <iostream>
#include <unordered_map>

#include "axom/inlet.hpp"

int main()
{
  auto lr = axom::inlet::cpp11_compat::make_unique<axom::inlet::LuaReader>();

  // Parse example input file
  lr->parseString("values = { [1] = 'start', [2] = 'stop', [3] = 'pause' }");

  axom::sidre::DataStore ds;

  // Initialize Inlet
  axom::inlet::Inlet inlet(std::move(lr), ds.getRoot());

  // Register the verifier, which will verify the array values
  auto& vals = inlet.getGlobalTable().addStringArray("values");
  vals.registerVerifier([](const axom::inlet::Table& table) -> bool {
    auto map = table.get<std::unordered_map<int, std::string>>();
    bool startFound = false;
    bool stopFound = false;
    for(auto p : map)
    {
      if(p.second == "start")
      {
        startFound = true;
        std::cout << "Found start at index " << p.first << std::endl;
      }
      else if(p.second == "stop")
      {
        stopFound = true;
        std::cout << "Found stop at index " << p.first << std::endl;
      }
    }
    return startFound && stopFound;
  });

  // We expect verfication to pass since values array has 3 elements
  inlet.verify() ? std::cout << "Verification passed\n"
                 : std::cout << "Verification failed\n";

  // Print contents of map
  std::unordered_map<int, std::string> map = inlet["values"];
  std::cout << "\nMap Contents:\n";
  for(auto p : map)
  {
    std::cout << p.first << " " << p.second << std::endl;
  }

  return 0;
}
