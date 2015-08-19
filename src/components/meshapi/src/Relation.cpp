/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


#include "Relation.hpp"


namespace asctoolkit {
namespace meshapi {

/**
 * \brief Definition of static instance of nullSet for all relations
 * \note Should this be a singleton or a global object?  Should the scope be public?
 */
  NullSet Relation::s_nullSet;

} // namespace meshapi
} // namespace asctoolkit
