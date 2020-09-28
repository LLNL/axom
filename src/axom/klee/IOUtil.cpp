// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/klee/IOUtil.hpp"

#include <stdexcept>

namespace axom { namespace klee { namespace internal {

std::vector<double> toDoubleVector(const conduit::Node &listNode,
        std::size_t expectedSize) {
    if (!listNode.dtype().is_number()) {
        std::ostringstream message;
        message << "Not an array of numbers: " << listNode.path();
        throw std::invalid_argument(message.str());
    }

    conduit::Node valueAsDoubleArray;
    listNode.to_double_array(valueAsDoubleArray);
    double *values = valueAsDoubleArray.as_double_ptr();
    std::vector<double> result(valueAsDoubleArray.dtype().number_of_elements());
    std::copy(values, values + result.size(), result.begin());

    if (result.size() != expectedSize) {
        std::ostringstream message;
        message << listNode.name() << " should be a list of " << expectedSize << " numbers";
        throw std::invalid_argument(message.str());
    }
    return result;
}

double toDouble(const conduit::Node &value) {
    auto &type = value.dtype();
    if (!type.is_number() || type.number_of_elements() != 1) {
        std::ostringstream message;
        message << value.name() << " should be a single number";
        throw std::invalid_argument(message.str());
    }
    return value.to_double();
}

}}}