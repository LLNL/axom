// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_NVTXRANGE_HPP_
#define AXOM_NVTXRANGE_HPP_

#include "axom/core/Macros.hpp"  // for axom macros

// C/C++ includes
#include <string> // for std::string

namespace axom
{

/*!
 * \brief Predefined set of NVTX colors to use with NVTXRange.
 */
enum class NVTXColor : uint32_t
{
  BLACK   = 0x00000000,
  GREEN   = 0x0000FF00,
  LIME    = 0x00BFFF00,
  RED     = 0x00FF0000,
  BLUE    = 0x000000FF,
  YELLOW  = 0x00FFFF00,
  CYAN    = 0x0000FFFF,
  MAGENTA = 0x00FF00FF,
  WHITE   = 0x00FFFFFF,
  ORANGE  = 0x00FFA500,
  PINK    = 0x00FF69B4
};

/*!
 * \brief Predefined set of user-specified categories to use with NVTXRange.
 */
enum class NVTXCategory : uint32_t
{
  ANY = 0,          /*!< wildcard used to annotate generic sections */
  PACKING,          /*!< used to annotate code that is doing packing */
  MEMTRANSFER,      /*!< used to annotate code that is doing memory transfer */
};

/*!
 * \brief Default NVTX color to use. Set to GREEN.
 */
constexpr NVTXColor DEFAULT_NVTX_COLOR = NVTXColor::GREEN;

/*!
 * \brief Default NVTX category to use. Set to ANY.
 */
constexpr NVTXCategory DEFAULT_NVTX_CATEGORY = NVTXCategory::ANY;

/*!
 * \class NVTXRange
 *
 * \brief NVTXRange is a simple utility class to annotate code.
 *
 *  The NVTXRange class is a simple utility class that can be used in
 *  conjunction with the NVIDIA Tools Extension library to allow developers
 *  to easily mark and annotate code in order to provide additional information
 *  to NVIDIA performance tools, such as, nvprof, nvvp and Nsight.
 *
 * \see https://docs.nvidia.com/cuda/profiler-users-guide/index.html#nvtx
 *
 * \note NVTXRange uses the RAII idiom, consequently the range is started
 *  when the NVTXRange object is instantiated and stopped when the object
 *  goes out of scope.
 *
 * \thanks Jason Burmark (burmark1@llnl.gov) for his original implementation
 *  that inspired the implementation of this class.
 *
 * Usage Example:
 * \code
 *
 *    // use scope to auto-start and stop a range
 *    { // begin scope resolution
 *    axom::NVTXRage range ("foo" );
 *    foor();
 *    } // end scope resoltuion
 *
 * \endcode
 *
 */
class NVTXRange
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  NVTXRange() = delete;

  /*!
   * \brief Creates an NVTXRage instance with the given name.
   *
   * \param [in] name the name to associate with the range
   * \param [in] color the color to associate with the range (optional)
   * \param [in] category the category to associate with the range (optional)
   *
   * \pre name.empty() == false
   */
  NVTXRange( const std::string& name,
             NVTXColor color = DEFAULT_NVTX_COLOR,
             NVTXCategory category= DEFAULT_NVTX_CATEGORY );

  /*!
   * \brief Destructor.
   */
  ~NVTXRange();

private:

  /*!
   * \brief Starts an NVTX range.
   * \note Called by the constructor.
   */
  void start();

  /*!
   * \brief Stops the NVTX range.
   * \note Called by the destructor.
   */
  void stop();

  std::string m_name;
  NVTXColor m_color;
  NVTXCategory m_category;
  bool m_active;

  DISABLE_COPY_AND_ASSIGNMENT(NVTXRange);
  DISABLE_MOVE_AND_ASSIGNMENT(NVTXRange);
};

} /* namespace axom */

#endif /* AXOM_NVTXRANGE_HPP_ */
