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
 * \brief Subsetting policies for the meshapi
 *
 * Subsetting policies encompass the type and availability of a set's parent
 *   A valid subset policy must support the following interface:
 *      [required]
 *      * indirection() : IntType  -- returns the value of the element after indirection
 *      * hasIndirection(): bool -- returns whether there is an indirection buffer
 *      [optional]
 *      * operator(): IntType -- alternate accessor for indirection
 */

#ifndef MESHAPI_POLICIES_SUBSET_H_
#define MESHAPI_POLICIES_SUBSET_H_


namespace asctoolkit {
namespace meshapi {
namespace policies {

    /**
     * \name OrderedSet_Subsetting_Policies
     * \brief A few default policies for the subsetting of an OrderedSet
     */

    /// \{

    struct NoSubset
    {
        static const NullSet s_nullSet;

        bool isSubset() const { return false; }
        const Set* parentSet() const { return &s_nullSet; }
    };

    struct VirtualParentSubset
    {
        static NullSet s_nullSet;

        VirtualParentSubset() : m_parentSet(&s_nullSet) {}

        bool isSubset() const { return *m_parentSet != s_nullSet; }
        const Set* parentSet() const { return m_parentSet; }
              Set*& parentSet()       { return m_parentSet; }

    private:
        Set* m_parentSet;
    };

    template<typename ParentSetType>
    struct ConcreteParentSubset
    {

        bool isSubset() const { return m_parentSet == ATK_NULLPTR; }
        const ParentSetType* parentSet() const { return m_parentSet; }
              ParentSetType*& parentSet()       { return m_parentSet; }

    private:
        ParentSetType* m_parentSet;
    };


    /// \}

} // end namespace policies
} // end namespace meshapi
} // end namespace asctoolkit

#endif // MESHAPI_POLICIES_SUBSET_H_
