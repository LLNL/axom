/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and
 * further review from Lawrence Livermore National Laboratory.
 */

/**
 * \file
 *
 * \brief Offset policies for the meshapi
 *
 * Offset policies are meant to represent the offset to the first element of meshapi ordered set
 *   A valid offset policy must support the following interface:
 *      [required]
 *      * DEFAULT_VALUE is a public static const IntType
 *      * offset() : IntType  -- returns the offset
 *      [optional]
 *      * operator(): IntType -- alternate accessor for the offset value
 */

#ifndef MESHAPI_POLICIES_OFFSET_H_
#define MESHAPI_POLICIES_OFFSET_H_


namespace asctoolkit {
namespace meshapi {
namespace policies {

    /**
     * \name OrderedSet_Offset_Policies
     * \brief A few default policies for the offset of an OrderedSet
     */

    /// \{

    /**
     * \brief A policy class for the offset in a set.  The offset can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeOffsetHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType();

        RuntimeOffsetHolder(IntType off = DEFAULT_VALUE) : m_off(off) {}

        inline IntType offset() const { return m_off;}
        inline IntType& offset()            { return m_off;}

        inline IntType  operator()() const { return offset();}
        inline IntType& operator()()       { return offset();}

        inline bool isValid(bool) const     { return true; }
    private:
        IntType m_off;
    };


    /**
     * \brief A policy class for a compile-time known set offset
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeOffsetHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeOffsetHolder(IntType val = DEFAULT_VALUE) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeOffsetHolder -- tried to initialize a compile time offset with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline IntType  offset()      const { return INT_VAL;}
        inline IntType operator()()  const { return offset();}
        inline bool isValid(bool) const     { return true; }
    };

    /**
     * \brief A policy class for when we have no offset
     */
    template<typename IntType>
    struct ZeroOffset {
        static const IntType DEFAULT_VALUE = IntType();

        ZeroOffset(IntType val = DEFAULT_VALUE) {
               SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                              , "MeshAPI::ZeroOffset policy -- tried to initialize a NoOffset policy with ("
                              << val <<", but should always be 0");
        }

        inline IntType offset()     const { return DEFAULT_VALUE;}
        inline IntType operator()() const { return offset();}
        inline bool isValid(bool) const     { return true; }
    };

    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_OFFSET_H_
