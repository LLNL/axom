// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#include "axom/sina.hpp"

int main(void) {
    // Define 3 different datums
    axom::sina::Datum myDatum{12.34};
    std::string value = "foobar";
    axom::sina::Datum myOtherDatum{value};
    std::vector<double> scalars = {1, 2, 20.0};
    axom::sina::Datum myArrayDatum{scalars};

    // Create a record to store the datum
    axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
    std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

    // Add the datum instances to the record
    myRecord->add("datum1", std::move(myDatum));
    myRecord->add("datum2", std::move(myOtherDatum));
    myRecord->add("datum3", std::move(myArrayDatum));

    // Query the datum
    auto &data = myRecord->getData();

    // Print the datum values
    std::cout << "datum1: " << data.at("datum1").getScalar() << std::endl;
    std::cout << "datum2: " << data.at("datum2").getValue() << std::endl;
    std::cout << "datum3: ";
    for (const auto& value : data.at("datum3").getScalarArray()) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}