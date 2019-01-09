/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * \file HandleMesh.cpp
 *
 * \brief Slam example where the ElementType of a set is a handle
 * that wraps a mesh index.
 */

#include "axom/slic.hpp"
#include "axom/slam.hpp"

#include <ostream>

namespace slam = axom::slam;
namespace slic = axom::slic;


namespace
{

/**
 * \class Handle
 * \brief A handle is a wrapper around a variable.
 *
 * This class is a stub for some mesh data structures
 * that use typed index handles rather than (untyped) indices
 * to refer to mesh elements.
 */
template<typename T>
struct Handle
{
  Handle() : mID(T()) {}
  Handle(T id) : mID(id) {}
  Handle(const Handle& h) : mID(h.mID) {}
  bool operator==(const Handle& h) { return mID == h.mID; }

  static Handle make_handle(T id) { return Handle(id); }

  void print(std::ostream& os) const
  {
    os << "Handle(" << mID <<")";
  }

  T mID;
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const Handle<T>& h)
{
  h.print(os);
  return os;
}

} // end anonymous namespace


int main(int, char**)
{
  slic::UnitTestLogger logger;

  using PosType = slam::PositionType;
  using HandleType = Handle<PosType>;

  using HandleSet = slam::VectorIndirectionSet<PosType, HandleType>;

  const int sz = 5;
  HandleSet::IndirectionBufferType vecHandle(sz);

  // Create a set of handles
  HandleSet hSet = HandleSet::SetBuilder()
                   .size( sz )
                   .data( &vecHandle );

  // Add handles with (somewhat) arbitrary IDs to the set
  for(auto i : hSet.positions() )
  {
    hSet[i] = HandleType::make_handle(2 * sz -i);
  }

  // Iterate over the set
  SLIC_INFO( "Iterating a set of Handle elements: " );
  for(auto i : hSet.positions() )
  {
    SLIC_INFO("  " << i << ": " << hSet[i] );
  }

  return 0;
}
