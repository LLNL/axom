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

        inline IndirectionResult indirection(PositionType pos) const { return static_cast<ElementType>(pos); }
        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        bool hasIndirection() const { return false;}
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct ArrayIndirection
    {
        typedef const ElementType& IndirectionResult;

        ArrayIndirection(ElementType* buf = ATK_NULLPTR) : m_arrBuf(buf) {}

        ElementType*& data() { return m_arrBuf;}

        inline IndirectionResult indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection()
                           , "MeshAPI::Set:ArrayIndirection -- Tried to dereference a null array in an array based indirection set.");
            return m_arrBuf[pos];
        }

        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        bool hasIndirection() const { return m_arrBuf == ATK_NULLPTR;}
    private:
        ElementType* m_arrBuf;
    };

    /**
     * \brief A policy class for sets with array-based indirection
     */
    template<typename PositionType, typename ElementType>
    struct STLVectorIndirection
    {
        typedef std::vector<ElementType> VectorType;
        typedef const ElementType& IndirectionResult;

        STLVectorIndirection(const VectorType* buf = ATK_NULLPTR) : m_vecBuf(buf) {}

        const VectorType*& data() { return m_vecBuf;}

        inline IndirectionResult indirection(PositionType pos) const
        {
            SLIC_ASSERT_MSG( hasIndirection(), "MeshAPI::Set:STLVectorIndirection -- Tried to dereference a null vector in a vector based indirection set.");
            //SLIC_ASSERT_MSG( pos < m_vecBuf->size(), "MeshAPI::Set:STLVectorIndirection -- Tried to access an out of bounds element at position "
            //        << pos << " in vector with only " << m_vecBuf->size() << " elements.");

            return (*m_vecBuf)[pos];
        }
        inline IndirectionResult operator()(PositionType pos)  const { return indirection(pos); }

        bool hasIndirection() const { return m_vecBuf != ATK_NULLPTR;}
    private:
        const VectorType* m_vecBuf;
    };

    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_INDIRECTION_H_
