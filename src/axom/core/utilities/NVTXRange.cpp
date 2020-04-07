// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/core/utilities/NVTXRange.hpp"

#include "axom/config.hpp"  // for axom compile-time definitions

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

NVTXRange::NVTXRange( const std::string& name,
                      NVTXColor color,
                      NVTXCategory category ) :
  m_name( name ),
  m_color( color ),
  m_category( category ),
  m_active( false )
{
  assert( !m_name.empty() );
  start();
  assert( m_active );
}

//------------------------------------------------------------------------------
NVTXRange::~NVTXRange()
{
  stop();
}

//------------------------------------------------------------------------------
void NVTXRange::start()
{
  assert( !m_active );
#ifdef AXOM_USE_CUDA

  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version       = NVTX_VERSION;
  eventAttrib.size          = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.category      = static_cast< uint32_t >( m_category );
  eventAttrib.colorType     = NVTX_COLOR_ARGB;
  eventAttrib.color         = static_cast< uint32_t>( m_color );
  eventAttrib.messageType   = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii = m_name.c_str();

  nvtxRangePushEx(&eventAttrib);

#endif
  m_active = true;
}

//------------------------------------------------------------------------------
void NVTXRange::stop()
{
  if ( m_active )
  {
#ifdef AXOM_USE_CUDA
    nvtxRangePop();
    m_active = false;
#endif
  }
}

} /* namespace axom */
