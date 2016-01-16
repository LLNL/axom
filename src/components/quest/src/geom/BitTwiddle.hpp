
#ifndef BIT_TWIDDLE_HXX_
#define BIT_TWIDDLE_HXX_

#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include <boost/type_traits.hpp>
# ifndef BOOST_NO_SFINAE
#  include <boost/utility/enable_if.hpp>
# endif

#include <limits>

/**
 * \file
 * \brief Some helper functions for efficient bitwise operations.
 *
 * Currently has some function to convert Points to/from Morton indices,
 * and a hash functor that uses these on Point types.
 *
 * \note The current implementation assumes that the Morton index is a 32-bit int
 *       (even though its type is std::size_t)
 *       This can easily be relaxed so we can handle deeper levels,
 *       but this is not yet implemented.  (It will require some type traits
 *       to ensure we have the right bit widths).
 *      Thus, for the time being, we assume that a 2D coordinate has at most 16 bits
 *      of information and a 3D coordinate has at most 10 bits of information.
 *      This can be extended to 32 and 21 bits, respectively if we want to fit into
 *      a 64 bit integer.
 */

namespace {

    /**
     * \brief Finds the index of the maximum bit set to one
     * An integer log2 function
     */
    template<typename CoordType>
    int maxSetBit(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        // Note, this is not a particularly fast implementation,
        // but is currently intended only for diagnostic purposes.

        const int maxBits = std::numeric_limits<CoordType>::digits;

        if(x == 0)
            return maxBits;

        int ret = 0;
        for(; x > 1 && ret < maxBits; ++ret, x >>=1) {}

        return ret;
    }

    // used for morton 2D and 4D
    unsigned int const B2[] = {
                static_cast< unsigned int >(0x5555555555555555),     // 0101'0101
                static_cast< unsigned int >(0x3333333333333333),     // 0011'0011
                static_cast< unsigned int >(0x0F0F0F0F0F0F0F0F),     // 0000'1111
                static_cast< unsigned int >(0x00FF00FF00FF00FF),     // x8
                static_cast< unsigned int >(0x0000FFFF0000FFFF),     // x16
                static_cast< unsigned int >(0x00000000FFFFFFFF), };  // x32

    // used for morton 2D and 4D (on short ints)
    unsigned short const sB2[] = {
                static_cast< unsigned short>(0x5555555555555555),     // 0101'0101
                static_cast< unsigned short>(0x3333333333333333),     // 0011'0011
                static_cast< unsigned short>(0x0F0F0F0F0F0F0F0F),     // 0000'1111
                static_cast< unsigned short>(0x00FF00FF00FF00FF),     // x8
                static_cast< unsigned short>(0x0000FFFF0000FFFF),     // x16
                static_cast< unsigned short>(0x00000000FFFFFFFF), };  // x32


    unsigned int const S2[] = {1, 2, 4, 8, 16, 32};

    // used for morton 3D
    // numbers from C. Ericson's Real Time Collision Detection book
    unsigned int const B3[] = {
                static_cast< unsigned int >(0x9249249249249249),     // 0010'0100'1001'0010'0100'1001
                static_cast< unsigned int >(0x30C30C30C30C30C3),     // 0000'1100'0011'0000'1100'0011
                static_cast< unsigned int >(0xF00F00F00F00F00F),     // 0000'0000'1111'0000'0000'1111
                static_cast< unsigned int >(0x00FF0000FF0000FF),     // 0000'0000'0000'0000'1111'1111
                static_cast< unsigned int >(0xFFFF00000000FFFF),     // x16
                static_cast< unsigned int >(0x00000000FFFFFFFF),     // x32
            };

    unsigned int const S3[] = {2, 4, 8, 16, 32, 64};


    template<typename CoordType>
    inline CoordType mortonize2(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        SLIC_ASSERT_MSG( x < (1<< 16), "Morton2D: Morton indexing in 2D currently only supports"
                         << " 16 bits per coordinate.  Attempted to index an integer (" << x
                         <<") with " <<  maxSetBit(x) << " bits.");

        x = (x | (x << S2[4])) & B2[4];
        x = (x | (x << S2[3])) & B2[3];
        x = (x | (x << S2[2])) & B2[2];
        x = (x | (x << S2[1])) & B2[1];
        return (x | (x << S2[0])) & B2[0];
    }

    inline short mortonize2(short x)
    {
        x = (x | (x << S2[3])) & sB2[3];
        x = (x | (x << S2[2])) & sB2[2];
        x = (x | (x << S2[1])) & sB2[1];
        return (x | (x << S2[0])) & sB2[0];
    }


    template<typename CoordType>
    inline CoordType demortonize2(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = (x | (x >> S2[0])) & B2[1];
        x = (x | (x >> S2[1])) & B2[2];
        x = (x | (x >> S2[2])) & B2[3];
        x = (x | (x >> S2[3])) & B2[4];
        return (x | (x >> S2[4])) & B2[5];

    }

    template<typename CoordType>
    inline CoordType demortonize3(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = (x | (x >> S3[0])) & B3[1];
        x = (x | (x >> S3[1])) & B3[2];
        x = (x | (x >> S3[2])) & B3[3];
        return (x | (x >> S3[3])) & B3[4];
    }



    template<typename CoordType>
    inline CoordType mortonize3(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);
        SLIC_ASSERT_MSG( x < (1<< 10), "Morton3D: Morton indexing in 3D currently only supports"
                         << " 10 bits per coordinate.  Attempted to index an integer (" << x
                         <<") with " <<  maxSetBit(x) << " bits.");

        x = (x | (x << S3[3])) & B3[3];
        x = (x | (x << S3[2])) & B3[2];
        x = (x | (x << S3[1])) & B3[1];
        return (x | (x << S3[0])) & B3[0];
    }
}

namespace quest
{
    typedef std::size_t MortonIndex;


    /**
     * \brief A function to convert a 2D point to a Morton index
     *
     * Morton indexing interleaves the bits of the point's coordinates
     * \param [in] x,y The coordinates of the point
     * \pre CoordType must be an integral type
     * \pre x and y must be positive with value less than \f$ 2^{16} \f$
     * \note These preconditions can easily be relaxed, if desired
     * (e.g. to 32 bits for 64 bit Morton indices)
     * \return The MortonIndex of the 2D point
     */
    template<typename CoordType>
    inline MortonIndex convertPointToMorton2D(CoordType x, CoordType y)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        return static_cast<MortonIndex>( mortonize2(x) | (mortonize2(y)<< 1) );
    }

    /**
     * \brief A function to convert a 3D point to a Morton index
     *
     *  Morton indexing interleaves the bits of the point's coordinates
     * \param [in] x,y,z The coordinates of the point
     * \pre CoordType must be an integral type
     * \pre x and y must be positive with value less than \f$ 2^{10} \f$
     * \note These preconditions can easily be relaxed, if desired
     * (e.g. to 21 bits for 64 bit Morton indices)
     * \return The MortonIndex of the 2D point
     */
    template<typename CoordType>
    inline MortonIndex convertPointToMorton3D(CoordType x, CoordType y, CoordType z)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        return static_cast<MortonIndex>( mortonize3(x) | (mortonize3(y) << 1) | (mortonize3(z) <<2 ) );
    }


    /**
     * \brief A function to convert a Morton index back to a 2D point
     *
     * \param [in] morton The MortonIndex of the desired point
     * \param [out] x,y The coordinates of the point
     *  Morton indexing interleaves the bits of the point's coordinates
     * \return The point's coordinates are in the x and y parameters
     */
    template<typename CoordType>
    inline void convertMortonToPoint2D(MortonIndex morton, CoordType &x, CoordType & y)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = demortonize2( morton      & B2[0]);
        y = demortonize2((morton >>1) & B2[0]);
    }


    /**
     * \brief A function to convert a Morton index back to a 3D point
     *
     * \param [in] morton The MortonIndex of the desired point
     * \param [out] x,y,z The coordinates of the point
     *  Morton indexing interleaves the bits of the point's coordinates
     * \return The point's coordinates are in the x, y and z parameters
     */
    template<typename CoordType>
    void convertMortonToPoint3D(MortonIndex morton, CoordType &x, CoordType & y, CoordType &z)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = demortonize3( morton      & B3[0]);
        y = demortonize3((morton >>1) & B3[0]);
        z = demortonize3((morton >>2) & B3[0]);
    }

    /**
     * \brief A function to convert a 2D point directly to a MortonIndex
     * \see convertMortonToPoint2D()
     *
     * \return The Morton index of the point
     */
    template<typename CoordType>
    inline MortonIndex convertPointToMorton(const Point<CoordType,2>& pt)
    {
        return convertPointToMorton2D(pt[0],pt[1]);
    }

    /**
     * \brief A function to convert a 3D point directly to a MortonIndex
     * \see convertMortonToPoint3D()
     *
     * \return The Morton index of the point
     */
    template<typename CoordType>
    inline MortonIndex convertPointToMorton(const Point<CoordType,3>& pt)
    {
        return convertPointToMorton3D(pt[0],pt[1], pt[2]);
    }

    /**
     * \brief A function to convert MortonIndex back to a 2D point
     * \see convertMortonToPoint2D()
     *
     * \return The demortonized 2D Point
     */
    template<typename CoordType>
    inline Point<CoordType,2> convertMortonToPoint2D(MortonIndex idx)
    {
        Point<CoordType,2> pt;
        convertMortonToPoint2D<CoordType>(idx, pt[0], pt[1]);
        return pt;
    }

    /**
     * \brief A function to convert MortonIndex back to a 3D point
     * \see convertMortonToPoint3D()
     *
     * \return The demortonized 3D Point
     */
    template<typename CoordType>
    inline Point<CoordType,3> convertMortonToPoint3D(MortonIndex idx)
    {
        Point<CoordType,3> pt;
        convertMortonToPoint3D(idx, pt[0], pt[1], pt[2]);
        return pt;
    }

    /**
     * \class
     * \brief A functor class for Mortonizing points.
     *
     * This can be used as a hashing function for Points in dimensions 1..4
     */
    template<typename CoordType>
    struct PointHash
    {
        /**
         * \brief Mortonizes a coordinate (viewed as a 1D point)
         * \note This is a no-op and is provided for genericity in point dimension
         * \param [in] coord The coordinate of the
         * \returns The morton index of the 1D point
         */
        std::size_t operator()(CoordType & coord) const
        {
            return coord;
        }

        /**
         * \brief Mortonizes a 1D point
         * \note This is a no-op and is provided for genericity in point dimension
         * \param [in] pt The 1D point
         * \returns The morton index of the point
         */
        std::size_t operator()(Point<CoordType,1> const& pt) const
        {
            return pt[0];
        }

        /**
         * \brief Mortonizes a 2D point
         * \param [in] pt The 2D point
         * \returns The morton index of the point
         */
        std::size_t operator()(Point<CoordType,2> const& pt) const
        {
            return convertPointToMorton2D<CoordType>(pt[0], pt[1]);
        }

        /**
         * \brief Mortonizes a 3D point
         * \param [in] pt The 3D point
         * \returns The morton index of the point
         */
        std::size_t operator()(Point<CoordType,3> const& pt) const
        {
            return convertPointToMorton3D<CoordType>(pt[0], pt[1], pt[2]);
        }

        /**
         * \brief Mortonizes a 4D point
         * \param [in] pt The 4D point
         * \returns The morton index of the point
         */
        std::size_t operator()(Point<CoordType,4> const& pt) const
        {
            MortonIndex ind1 = convertPointToMorton2D<CoordType>(pt[0], pt[2]);
            MortonIndex ind2 = convertPointToMorton2D<CoordType>(pt[1], pt[3]);
            return convertPointToMorton2D<MortonIndex>(ind1,ind2);
        }
    };

} // end namespace quest

#endif  // BIT_TWIDDLE_HXX_
