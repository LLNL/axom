
#ifndef BIT_TWIDDLE_HXX_
#define BIT_TWIDDLE_HXX_

#include "quest/Point.hpp"
#include "quest/Vector.hpp"

#include <boost/type_traits.hpp>
# ifndef BOOST_NO_SFINAE
#  include <boost/utility/enable_if.hpp>
# endif

namespace {

    // used for morton 2D and 4D
    unsigned int const B2[] = {
                static_cast< unsigned int >(0x5555555555555555),     // 0101
                static_cast< unsigned int >(0x3333333333333333),     // 0011
                static_cast< unsigned int >(0x0F0F0F0F0F0F0F0F),     // 00001111
                static_cast< unsigned int >(0x00FF00FF00FF00FF),     // x8
                static_cast< unsigned int >(0x0000FFFF0000FFFF),     // x16
                static_cast< unsigned int >(0x00000000FFFFFFFF), };  // x32

    // used for morton 2D and 4D (on short ints)
    unsigned short const sB2[] = {
                static_cast< unsigned short>(0x5555555555555555),     // 0101
                static_cast< unsigned short>(0x3333333333333333),     // 0011
                static_cast< unsigned short>(0x0F0F0F0F0F0F0F0F),     // 00001111
                static_cast< unsigned short>(0x00FF00FF00FF00FF),     // x8
                static_cast< unsigned short>(0x0000FFFF0000FFFF),     // x16
                static_cast< unsigned short>(0x00000000FFFFFFFF), };  // x32


    unsigned int const S2[] = {1, 2, 4, 8, 16};

    // used for morton 3D
    // from C. Ericson's Real Time Collision Detection book
    //tested only in encoding direction (not decoding)
    unsigned int const B3[] = {
                static_cast< unsigned int >(0x9249249249249249),     // 001001001001001001001001
                static_cast< unsigned int >(0x30C30C30C30C30C3),     // 000011000011000011000011
                static_cast< unsigned int >(0xF00F00F00F00F00F),     // 000000001111000000001111
                static_cast< unsigned int >(0x00FF0000FF0000FF),     // 000000000000000011111111
            };

    unsigned int const S3[] = {2, 4, 8, 16};


    template<typename CoordType>
    inline CoordType mortonize2(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

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

        x = (x | (x >> S2[0])) & B2[0];
        x = (x | (x >> S2[1])) & B2[1];
        x = (x | (x >> S2[2])) & B2[2];
        x = (x | (x >> S2[3])) & B2[3];
        return (x | (x >> S2[4])) & B2[4];

    }

    template<typename CoordType>
    inline CoordType demortonize3(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = (x | (x >> S3[0])) & B3[0];
        x = (x | (x >> S3[1])) & B3[1];
        x = (x | (x >> S3[2])) & B3[2];
        return (x | (x >> S3[3])) & B3[3];
    }



    template<typename CoordType>
    inline CoordType mortonize3(CoordType x)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = (x | (x << S3[3])) & B3[3];
        x = (x | (x << S3[2])) & B3[2];
        x = (x | (x << S3[1])) & B3[1];
        return (x | (x << S3[0])) & B3[0];
    }
}

namespace quest
{
    typedef std::size_t MortonIndex;


    template<typename CoordType>
    inline MortonIndex convertPointToMorton2D(CoordType x, CoordType y)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        return static_cast<MortonIndex>( mortonize2(x) | (mortonize2(y)<< 1) );
    }

    template<typename CoordType>
    inline MortonIndex convertPointToMorton3D(CoordType x, CoordType y, CoordType z)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        return static_cast<MortonIndex>( mortonize3(x) | (mortonize3(y) << 1) | (mortonize3(z) <<2 ) );
    }

    template<typename CoordType>
    inline void convertMortonToPoint2D(MortonIndex morton, CoordType &x, CoordType & y)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = demortonize2( morton      & B2[0]);
        y = demortonize2((morton >>1) & B2[0]);
    }


    /// WARNING -- Untested, needs to be converted to parallel bitwise ops
    template<typename CoordType>
    void convertMortonToPoint3D(MortonIndex morton, CoordType &x, CoordType & y, CoordType &z)
    {
        BOOST_STATIC_ASSERT(boost::is_integral<CoordType>::value);

        x = demortonize3( morton      & B3[0]);
        y = demortonize3((morton >>1) & B3[0]);
        z = demortonize3((morton >>2) & B3[0]);
    }


    template<typename CoordType>
    inline MortonIndex convertPointToMorton(const Point<CoordType,2>& pt)
    {
        return convertPointToMorton2D(pt[0],pt[1]);
    }

    template<typename CoordType>
    inline MortonIndex convertPointToMorton(const Point<CoordType,3>& pt)
    {
        return convertPointToMorton3D(pt[0],pt[1], pt[2]);
    }

    template<typename CoordType>
    inline Point<CoordType,2> convertMortonToPoint(const MortonIndex idx)
    {
        Point<CoordType,2> pt;
        convertMortonToPoint2D(idx, &pt[0], &pt[1]);
        return pt;
    }

    template<typename CoordType>
    inline Point<CoordType,3> convertMortonToPoint(const MortonIndex idx)
    {
        Point<CoordType,3> pt;
        convertMortonToPoint2D(idx, &pt[0], &pt[1], &pt[2]);
        return pt;
    }

    template<typename CoordType>
    struct PointHash
    {
        std::size_t operator()(CoordType & coord) const
        {
            return coord;
        }

        std::size_t operator()(Point<CoordType,1> const& pt) const
        {
            return pt[0];
        }

        std::size_t operator()(Point<CoordType,2> const& pt) const
        {
            return convertPointToMorton2D(pt[0], pt[1]);
        }

        std::size_t operator()(Point<CoordType,3> const& pt) const
        {
            return convertPointToMorton3D(pt[0], pt[1], pt[2]);
        }

        std::size_t operator()(Point<CoordType,4> const& pt) const
        {
            MortonIndex ind1 = convertPointToMorton2D(pt[0], pt[2]);
            MortonIndex ind2 = convertPointToMorton2D(pt[1], pt[3]);
            return convertPointToMorton2D(ind1,ind2);
        }
    };

} // end namespace quest

#endif  // BIT_TWIDDLE_HXX_
