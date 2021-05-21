// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/nvtx/Range.hpp"

#include "axom/config.hpp"  // for axom compile-time definitions
#include "axom/core/utilities/nvtx/interface.hpp"

// C/C++ includes
#include <cassert>

// CUDA NVTX includes
#ifdef AXOM_USE_CUDA
  #include <cuda.h>
  #include <nvToolsExt.h>
  #include <nvToolsExtCuda.h>
#endif

namespace axom
{
namespace nvtx
{
Range::Range(const std::string& name) : m_name(name), m_active(false)
{
  assert(!m_name.empty());
  start();
  assert(m_active);
}

//------------------------------------------------------------------------------
Range::~Range() { stop(); }

//------------------------------------------------------------------------------
void Range::start()
{
  assert(!m_active);

#ifdef AXOM_USE_CUDA
  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version = NVTX_VERSION;
  eventAttrib.size = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.category = nvtx::get_category();
  eventAttrib.colorType = NVTX_COLOR_ARGB;
  eventAttrib.color = static_cast<uint32_t>(nvtx::get_color());
  eventAttrib.messageType = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = m_name.c_str();

  nvtxRangePushEx(&eventAttrib);
#endif

  m_active = true;
}

//------------------------------------------------------------------------------
void Range::stop()
{
  if(m_active)
  {
#ifdef AXOM_USE_CUDA
    nvtxRangePop();
#endif
    m_active = false;
  }
}

} /* namespace nvtx */

} /* namespace axom */
