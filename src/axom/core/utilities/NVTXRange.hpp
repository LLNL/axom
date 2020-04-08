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

// Forward Declarations
class NVTXRange;

/// \name AXOM NVTX Macros
///@{

/*!
 * \def AXOM_NVTX_SECTION
 *
 * \brief The AXOM_NVTX_SECTION macro is used to annotate sections of code
 *
 * \note In contrast to the AXOM_NVTX_FUNCTION macro, the AXOM_NVTX_SECTION
 *   macro is used to annotate sections of code, at a much finer granularity,
 *   within a given function.
 *
 * \warning Variables declared within a given AXOM_NVTX_SECTION are only defined
 *  within the scope of the AXOM_NVTX_SECTION.
 *
 * \warning An AXOM_NVTX_SECTION cannot be called in a nested fashion, i.e.,
 *  within another AXOM_NVTX_SECTION
 *
 * \note You may have multiple AXOM_NVTX_SECTION defined within a function and
 *  this macro can be used in conjunction with the AXOM_NVTX_FUNCTION macro.
 *
 * \Usage Example:
 * \code
 *
 *   void foo( )
 *   {
 *     AXOM_NVTX_FUNCTION( "foo"" );
 *
 *     // STEP 0: Run kernel A
 *     AXOM_NVTX_SECTION( "kernelA",
 *
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType i)
 *        {
 *          ..
 *        } );
 *
 *     ); // END NVTX SECTION for kernel A
 *
 *     // STEP 1: Run kernel B
 *     AXOM_NVTX_SECTION( "kernelB",
 *
 *        axom::for_all( 0, N, AXOM_LAMBDA(axom::IndexType i)
 *        {
 *          ...
 *        } );
 *
 *     ); // END NVTX SECTION for kernel B
 *
 *   }
 * \endcode
 *
 */
#if defined(AXOM_ENABLE_ANNOTATIONS) && defined(AXOM_USE_CUDA)
#define AXOM_NVTX_SECTION( __name__, ... )                                    \
  do {                                                                        \
    axom::NVTXRange r(__name__);                                              \
    __VA_ARGS__                                                               \
  } while( false )
#else
#define AXOM_NVTX_SECTION( __name__, ... )                                    \
  do {                                                                        \
    __VA_ARGS__                                                               \
  } while( false )
#endif

/*!
 * \def AXOM_NVTX_FUNCTION( name )
 *
 * \brief The AXOM_NVTX_FUNCTION macro is used to annotate a function.
 * \param [in] name a user-supplied name that will be given to the range.
 *
 * \note Typically, the AXOM_NVTX_FUNCTION macro is placed in the beginning of
 *  the function to annotate.
 *
 * \warning The AXOM_NVTX_FUNCTION can be called once within a (function) scope.
 *
 * Usage Example:
 * \code
 *   void foo( )
 *   {
 *     AXOM_NVTX_FUNCTION( "foo" );
 *     ...
 *   }
 * \endcode
 *
 */
#if defined(AXOM_ENABLE_ANNOTATIONS) && defined(AXOM_USE_CUDA)
#define AXOM_NVTX_FUNCTION( __name__ ) axom::NVTXRange __func_range(__name__)
#else
#define AXOM_NVTX_FUNCTION( __name__ )
#endif

///@}

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
