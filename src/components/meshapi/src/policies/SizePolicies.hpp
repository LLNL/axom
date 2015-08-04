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
 * \brief Size policies for the meshapi
 *
 * Size policies are meant to represent the size of a meshapi entity (e.g. the size of a set).
 *   A valid size policy must support the following interface:
 *      [required]
 *      * DEFAULT_VALUE is a public static const IntType
 *      * size() : IntType  -- returns the underlying integer size
 *      * empty() : bool -- returns whether the size is zero
 *      [optional]
 *      * operator(): IntType -- alternate accessor for the size value
 */

#ifndef MESHAPI_POLICIES_SIZE_H_
#define MESHAPI_POLICIES_SIZE_H_


namespace asctoolkit {
namespace meshapi {
namespace policies {

    /**
     * \name OrderedSet_Size_Policies
     * \brief A few default policies for the size of an OrderedSet
     */

    /// \{

    /**
     * \brief A policy class for the size of a set whose value can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeSizeHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType();

        RuntimeSizeHolder(IntType sz = DEFAULT_VALUE) : m_sz(sz) {}

        inline const IntType  size()       const { return m_sz;}
        inline       IntType& size()             { return m_sz;}

        inline const IntType  operator()() const { return size();}
        inline       IntType& operator()()       { return size();}

        inline bool empty() const       { return m_sz == IntType(); }
    private:
        IntType m_sz;
    };

    /**
     * \brief A policy class for a compile-time known set size
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeSizeHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeSizeHolder(IntType val = DEFAULT_VALUE) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeSizeHolder -- tried to initialize a compile time size policy with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline const IntType size()       const { return INT_VAL;}
        inline const IntType operator()() const { return size();}

        inline IntType empty() const { return INT_VAL == IntType();}
    };

    /**
      * \brief A policy class for an empty set (no size)
      */
     template<typename IntType>
     struct ZeroSize {
         static const IntType DEFAULT_VALUE = IntType();

         ZeroSize(IntType val = DEFAULT_VALUE) {
                SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                               , "MeshAPI::ZeroSize policy-- tried to initialize a NoSize set with value with value ("
                              << val <<" ) but should always be zero.");
         }

         inline const IntType size()        const { return DEFAULT_VALUE;}
         inline const IntType operator()()  const { return size();}
         inline IntType empty() const { return true;}
     };

    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_SIZE_H_
