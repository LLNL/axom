// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Create the record
  axom::sina::ID myID {"my_record", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> myRecord {
    new axom::sina::Record {myID, "my_type"}};

  // Create the datum with an array of strings
  std::vector<std::string> myTags {"input"};
  axom::sina::Datum myDatum {myTags};

  // Add the datum to the record
  myRecord->add("my_scalar", std::move(myDatum));
  std::cout << myRecord->toNode().to_json() << std::endl;
}