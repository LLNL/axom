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
 * \brief Indirection policies for the meshapi
 *
 * Indirection policies encompass the underlying storage for indirection buffers for a meshapi set, relation or map.
 *
 * A valid indirection policy must support the following interface:
 *      * typedef IndirectionResult -- the type of the result of an indirection (const/nonconst and ref/nonref)
 *      * indirection() : IntType  -- returns the value of the element after indirection
 *      * hasIndirection(): bool -- returns whether there is an indirection buffer
 *      * [optional] operator(): IntType -- alternate accessor for indirection
 *      * [optional] data() : ElementType* -- allows direct access to the underlying buffer (when this exists)
 */

#ifndef MESHAPI_POLICIES_INDIRECTION_H_
#define MESHAPI_POLICIES_INDIRECTION_H_

#include "slic/slic.hpp"




namespace asctoolkit {
namespace meshapi {
namespace policies {

    /**
     * \name OrderedSet_Indirection_Policies
     * \brief A few default policies for the indirection of an OrderedSet
     */

    /// \{

    /**
     * \brief A policy class for sets with no indirection
     */
    template<typename PositionType, typename ElementType>
    struct NoIndirection
    {
        typedef const ElementType IndirectionResult;
        typedef struct{} IndirectionBufferType;

        NoIndirection() {}
        NoIndirection(IndirectionBufferType*) {}     // This empty .ctor exists to satisfy IndirectionPolicy interface

        inline IndirectionResult indirection(PositionType pos) const { return static_cast<ElementType>(pos); }
        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        IndirectionBufferType* data() { return ATK_NULLPTR;}

        bool hasIndirection() const { return false;}
        inline bool isValid(PositionType, PositionType, PositionType, bool ) const     { return true; }
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct ArrayIndirection
    {
        typedef const ElementType& IndirectionResult;
        typedef ElementType IndirectionBufferType;

        ArrayIndirection(IndirectionBufferType* buf = ATK_NULLPTR) : m_arrBuf(buf) {}

        IndirectionBufferType*& data() { return m_arrBuf;}

        inline IndirectionResult indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection()
                           , "MeshAPI::Set:ArrayIndirection -- Tried to dereference a null array in an array based indirection set.");
            return m_arrBuf[pos];
        }

        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        bool hasIndirection() const { return m_arrBuf != ATK_NULLPTR;}

        inline bool isValid(PositionType size, PositionType, PositionType, bool verboseOutput = false) const
        {
            // set of zero size is always valid
            if(size == 0)
                return true;

            bool bValid = hasIndirection();
            SLIC_CHECK_MSG(verboseOutput && !bValid
                         , "Array-based indirection set with non-zero size (size=" << size
                           << ") requires valid data buffer, but buffer pointer was null.");

            // Since array-based indirection sets don't encode a buffer size, that is all we can check

            return bValid;
        }

    private:
        IndirectionBufferType* m_arrBuf;
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct STLVectorIndirection
    {
        typedef std::vector<ElementType> VectorType;
        typedef const ElementType& IndirectionResult;
        typedef const VectorType IndirectionBufferType;


        STLVectorIndirection(IndirectionBufferType* buf = ATK_NULLPTR) : m_vecBuf(buf) {}

        IndirectionBufferType*& data() { return m_vecBuf;}

        inline IndirectionResult indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection(), "MeshAPI::Set:STLVectorIndirection -- Tried to dereference a null vector in a vector based indirection set.");
            //SLIC_ASSERT_MSG( pos < m_vecBuf->size(), "MeshAPI::Set:STLVectorIndirection -- Tried to access an out of bounds element at position "
            //        << pos << " in vector with only " << m_vecBuf->size() << " elements.");

            return (*m_vecBuf)[pos];
        }
        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        bool hasIndirection() const { return m_vecBuf != ATK_NULLPTR;}

        inline bool isValid(PositionType size, PositionType offset, PositionType stride, bool verboseOutput = false) const
        {
            // If set has zero size, we are always valid (even if indirection buffer is null)
            if(size == 0)
                return true;

            // Otherwise, check whether the set has elements, but the array ptr is null
            bool bValid = hasIndirection();
            SLIC_CHECK_MSG(!verboseOutput || bValid
                          , "Vector-based indirection set with non-zero size (size=" << size
                              << ") requires valid data buffer, but buffer pointer was null.");
            if(!bValid)
                return false;

            // Finally, check that the underlying vector has sufficient storage for all set elements
            // Note that it is valid for the data buffer to have more space than the set's positions
            PositionType firstElt = offset;
            PositionType lastElt = (size-1)*stride + offset;
            PositionType vecSize = m_vecBuf->size();

            bValid = (firstElt < vecSize) && (lastElt < vecSize);
            SLIC_CHECK_MSG(!verboseOutput || bValid
                          , "Data buffer in vector-based IndirectionSet must be large enough to hold all elements of the set. "
                          << "Underlying buffer size is " << vecSize << ", and set's range is from "
                          <<"positions " << firstElt << " to " << lastElt <<".");

            return bValid;
        }

    private:
        IndirectionBufferType* m_vecBuf;
    };

    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_INDIRECTION_H_
