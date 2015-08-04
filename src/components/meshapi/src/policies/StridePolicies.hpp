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
 * \brief Stride policies for the meshapi
 *
 * Stride policies are meant to represent the fixed distance between consecutive elements of an OrderedSet
 *   A valid stride policy must support the following interface:
 *      [required]
 *      * DEFAULT_VALUE is a public static const IntType
 *      * stride() : IntType  -- returns the stride
 *      [optional]
 *      * operator(): IntType -- alternate accessor for the stride value
 */

#ifndef MESHAPI_POLICIES_STRIDE_H_
#define MESHAPI_POLICIES_STRIDE_H_


namespace asctoolkit {
namespace meshapi {
namespace policies {

    /**
     * \name OrderedSet_Stride_Policies
     * \brief A few default policies for the stride of an OrderedSet
     */

    /// \{


    /**
     * \brief A policy class for the stride in a set.  When using this class, the stride can be set at runtime.
     */
    template<typename IntType>
    struct RuntimeStrideHolder {
    public:
        static const IntType DEFAULT_VALUE = IntType(1);

        RuntimeStrideHolder(IntType stride = DEFAULT_VALUE) : m_stride(stride) {}

        inline const IntType  stride()      const { return m_stride;}
        inline       IntType& stride()            { return m_stride;}

        inline const IntType  operator()()  const { return stride();}
        inline       IntType& operator()()        { return stride();}


        //inline bool hasStride() const       { return m_stride != IntType(); }
    private:
        IntType m_stride;
    };


    /**
     * \brief A policy class for a compile-time known stride
     */
    template<typename IntType, IntType INT_VAL>
    struct CompileTimeStrideHolder {
        static const IntType DEFAULT_VALUE = INT_VAL;

        CompileTimeStrideHolder(IntType val = DEFAULT_VALUE) {
               SLIC_ASSERT_MSG( val == INT_VAL
                              , "MeshAPI::CompileTimeStrideHolder -- tried to initialize a compile time stride with value ("
                              << val <<" ) that differs from the template parameter of " << INT_VAL <<".");
        }

        inline const IntType stride()      const { return INT_VAL;}
        inline const IntType operator()()  const { return stride();}
    };

    /**
     * \brief A policy class for a set with stride one (i.e. the default stride)
     */
    template<typename IntType>
    struct StrideOne {
        static const IntType DEFAULT_VALUE = IntType(1);

        /**
         * This constructor only exists to allow the derived class to not have to specialize for when the stride is known at compile time
         */
        StrideOne(IntType val = DEFAULT_VALUE) {
               SLIC_ASSERT_MSG( val == DEFAULT_VALUE
                              , "MeshAPI::StrideOne policy -- tried to initialize a stride-one StridePolicy with value ("
                              << val <<"), but should always be 1.");
        }

        inline const IntType stride()     const { return DEFAULT_VALUE;}
        inline const IntType operator()() const { return stride();}
    };


    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_STRIDE_H_
