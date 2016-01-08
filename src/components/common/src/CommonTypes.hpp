/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/*!
 *  \file CommonTypes.hpp
 *
 *  \brief File containing common types used in toolkit components.
 *
 */

#ifndef COMMONTYPES_HPP_
#define COMMONTYPES_HPP_

#ifndef USE_CXX11
   #include <cstddef>       // brings in NULL
#endif


namespace asctoolkit
{

namespace common
{

#ifdef USE_CXX11
#define ATK_NULLPTR nullptr
#else
#define ATK_NULLPTR NULL
#endif


} /* end namespace common */
} /* end namespace asctoolkit */


#endif /* COMMONTYPES_HPP_ */
