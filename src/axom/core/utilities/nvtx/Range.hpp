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
namespace nvtx
{

/*!
 * \brief Predefined set of NVTX colors to use with NVTXRange.
 */
enum class Color : uint32_t
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
 * \brief Default NVTX color to use. Set to GREEN.
 */
constexpr Color DEFAULT_COLOR = Color::GREEN;

/*!
 * \brief Wildcard used for category
 */
constexpr uint32_t ANY_CATEGORY = 0;

/*!
 * \brief Default NVTX category to use. Set to ANY.
 */
constexpr uint32_t DEFAULT_CATEGORY = ANY_CATEGORY;


/*!
 * \class Range
 *
 * \brief Range is a simple utility class to annotate code.
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
 * \remark Thanks to Jason Burmark (burmark1@llnl.gov) for his original 
 *  implementation that inspired the implementation of this class.
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
class Range
{
public:

  /*!
   * \brief Default constructor. Disabled.
   */
  Range() = delete;

  /*!
   * \brief Creates an NVTXRage instance with the given name.
   *
   * \param [in] name the name to associate with the range
   * \param [in] color the color to associate with the range (optional)
   * \param [in] category the category to associate with the range (optional)
   *
   * \pre name.empty() == false
   */
  Range( const std::string& name,
         Color color = DEFAULT_COLOR,
         uint32_t category= DEFAULT_CATEGORY );

  /*!
   * \brief Destructor.
   */
  ~Range();

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
  Color m_color;
  uint32_t m_category;
  bool m_active;

  DISABLE_COPY_AND_ASSIGNMENT(Range);
  DISABLE_MOVE_AND_ASSIGNMENT(Range);
};

} /* namespace nvtx */

} /* namespace axom */

#endif /* AXOM_NVTXRANGE_HPP_ */
