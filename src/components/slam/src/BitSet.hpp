/*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
*
* Produced at the Lawrence Livermore National Laboratory
*
* LLNL-CODE-741217
*
* All rights reserved.
*
* This file is part of Axom.
*
* For details about use and distribution, please read axom/LICENSE.
*
*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/**
* \file BitSet.hpp
*
* \brief Contains a BitSet class for manipulating ordered sequences of bits.
*/

#ifndef SLAM_BITSET_H_
#define SLAM_BITSET_H_

#include "axom/Types.hpp"
#include "axom_utils/Utilities.hpp"
#include "slic/slic.hpp"

#include "slam/BitTwiddle.hpp"

#include <vector>

namespace axom
{
namespace slam
{

// Forward declare the BitSet class and some operator functions

class BitSet;

/**
 * \brief Union operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The union of the two input bitsets, where the i^th bit is set
 * it is set in either lhs or rhs (e.g. b[i] = lhs[i] || rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator|(const BitSet & lhs, const BitSet & rhs);

/**
 * \brief Intersection operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The intersection of the two input bitsets, where the i^th bit is set
 * it is set in both lhs and rhs (e.g. b[i] = lhs[i] && rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator&(const BitSet & lhs, const BitSet & rhs);

/**
 * \brief Exclusive or (xor) operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The xor of the two input bitsets, where the i^th bit is set
 * it is set in either lhs or rhs but not both (e.g. b[i] = lhs[i] ^ rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator^(const BitSet & lhs, const BitSet & rhs);

/**
 * \brief Set difference operator for two bit sets.
 *
 *
 * \param lhs The first bitset
 * \param rhs The second bitset
 * \return The set difference of the two input bitsets, where the i^th bit is
 * set it is set in lhs and not rhs (e.g. b[i] = lhs[i] && ! rhs[i]).
 *
 * \pre lhs.size() == rhs.size()
 */
BitSet operator-(const BitSet & lhs, const BitSet & rhs);

/**
* \class BitSet
* \brief A bitset class
*
* This class supports bitwise manipulation operations (e.g. set intersection,
* union and difference) on an ordered set of bits. The class has a similar
* interface to std::bitset and boost::dynamic_bitset, but with the following
* differences:
*   - The size of the bitset is supplied at runtime in the class constructor.
*   However, BitSet does not currently support changing the size after construction.
*   - We do not support the random access operation ( operator[](index) ).
*   The value of individual bits can be checked using the test function
*   (e.g. bset.test(i) )
*   - There is no support for directly initializing the bits in the bitset
*   (e.g. via strings).
*
*   The individual bits in the bitset are packed into a contiguous array of Words
*   (an unsigned integer type). Many bitset operations such as count() and flip()
*   are performed at the granularity of a word.
*
*   BitSet uses the same interface as boost::dynamic_bitset to enumerate the
*   set bits. Specifically, find_first() returns the index of the first set bit,
*   and find_next(idx) returns the next bit that is set after bit index idx.
*   BitSet::npos is used as a sentinal to indicate that there are no more set bits.
*/
class BitSet
{
public:
    typedef int Index;
    typedef axom::common::uint64 Word;

    // Use vector for initial implementation -- TODO: update using a policy
    typedef std::vector<Word> ArrayType;

    static const Index npos;

private:

    enum
    {
        BITS_PER_WORD = internal::BitTraits<Word>::BITS_PER_WORD,
        LG_BITS_PER_WORD = internal::BitTraits<Word>::LG_BITS_PER_WORD
    };

public:

    /**
     * \brief BitSet class constructor
     *
     * \param numBits The number of bits in the bitset.
     * \pre numBits must be non-negative
     * \post bset.size() == numBits
     * \post All bits will be off
     */
    explicit BitSet(int numBits = 0)
    {
        SLIC_ASSERT_MSG(numBits >= 0,
            "slam::BitSet must be initialized with a non-zero number of bits");

        m_numBits = axom::utilities::max(numBits, 0);
        m_numWords = (m_numBits == 0)
            ? 1
            : 1 + (m_numBits - 1) / BITS_PER_WORD;

        m_data = ArrayType(m_numWords);
    }

    /** Copy constructor for BitSet class */
    BitSet(const BitSet& other)
        : m_data(other.m_data),
          m_numBits(other.m_numBits),
          m_numWords(other.m_numWords)
    {}

    /** Assignment operator for BitSet class */
    BitSet& operator=(const BitSet& other)
    {
        if (this != &other)
        {
            m_data = other.m_data;
            m_numBits = other.m_numBits;
            m_numWords = other.m_numWords;
        }
        return *this;
    }
    
    /**
     * BitSet union-assignment operator
     *
     * \param other The other bitset
     * \pre other.size() == size()
     * \return A reference to the modified bitset
     *
     * Set union of the current instance and other. The i^th
     * bit will be set if it is set in this bitset or in \a other
     */
    BitSet& operator|=(const BitSet& other);

    /**
     * BitSet intersection-assignment operator
     *
     * \param other The other bitset
     * \pre other.size() == size()
     * \return A reference to the modified bitset
     *
     * Set intersection of the current instance and other. The i^th
     * bit will be set if it is set in this bitset and in \a other
     */
    BitSet& operator&=(const BitSet& other);

    /**
     * BitSet exclusive-or-assignment operator
     *
     * \param other The other bitset
     * \pre other.size() == size()
     * \return A reference to the modified bitset
     *
     * Set exclusive-or of the current instance and other. The i^th
     * bit will be set if it is set in this bitset or in \a other,
     * but not both.
     */
    BitSet& operator^=(const BitSet& other);

    /**
     * BitSet difference-assignment operator
     *
     * \param other The other bitset
     * \pre other.size() == size()
     * \return A reference to the modified bitset
     *
     * Set difference of the current instance and other. The i^th
     * bit will be set if it is set in this bitset and not in \a other
     */
    BitSet& operator-=(const BitSet& other);

    /** Equality operator for two bit sets */
    bool operator==(const BitSet & other) const;

    /** Inequality operator for two bit sets */
    bool operator!=(const BitSet & other) const
    {
        return !(operator==(other));
    }

public:
    /// \name Bitset iteration interface
    /// @{

    /**
     * Finds the index of the first bit that is set in the bitset
     *
     * \return The index of the first set bit,
     * or BitSet::npos if no bits are set
     */
    Index find_first() const;

    /**
     * Finds the index of the next bit that is set in the bitset after \a idx
     *
     * \param idx The starting index
     * \return The index of the first set bit after index \a idx
     * or BitSet::npos if none can be found
     * \note Will also return BitSet::npos if \a idx is BitSet::npos
     */
    Index find_next(Index idx) const;

    /// @}

public:
    /// \name Operations that affect all bits in the bitset
    /// @{

    /** Returns the cardinality of the bitset */
    int size() const { return m_numBits; }

    /** Returns the number of bits that are set */
    int count() const;

    /** Clears all bits in the bitset */
    void clear();

    /** Sets all bits in the bitset */
    void set();

    /** Toggles all bits in the bitset */
    void flip();

    /**
     * Checks if the bitset instance is valid
     *
     * \return True if the bitset is valid, false otherwise
     *
     * A bitset is valid if it has sufficient storage for size() bits.
     * If we have storage for more than size() bits, none of these additional
     * bits are set.
     */
    bool isValid() const;

    /// @}

public:
    /// \name Operations that affect a single bits in the bitset
    /// @{

    /**
     * Clears bit at index \a idx
     *
     * \pre \a idx must be between 0 and bitset.size()
     */
    void clear(Index idx) { getWord(idx) &= ~mask(idx); }

    /**
     * Sets bit at index \a idx
     *
     * \pre \a idx must be between 0 and bitset.size()
     */
    void set(Index idx)   { getWord(idx) |= mask(idx);  }

    /**
     * Toggles bit at index \a idx
     *
     * \pre \a idx must be between 0 and bitset.size()
     */
    void flip(Index idx)  { getWord(idx) ^= mask(idx);  }

    /**
     * Tests the bit at index \a idx
     *
     * \return True if bit \idx is set, false otherwise
     * \pre \a idx must be between 0 and bitset.size()
     */
    bool test(Index idx)  const
    {
        return (getWord(idx) & mask(idx)) != Word(0);
    }

    /// @}
private:
    /**
     * Gets the index of the word containing bit at index \a idx
     *
     * \param idx The index of the word containing the desired bit
     * \param checkIndexValid Option to enable bounds checking to ensure
     * that \a idx is within range [0, size() )
     */
    Word& getWord(Index idx, bool checkIndexValid = true)
    {
        if(checkIndexValid)
            checkValidIndex(idx);

        const Index wIdx = idx / BITS_PER_WORD;
        return m_data[wIdx];
    }

    /**
     * Const implementation of getWord()
     * \sa getWord()
     */
    const Word& getWord(Index idx, bool checkIndexValid = true) const
    {
        if (checkIndexValid)
            checkValidIndex(idx);

        const Index wIdx = idx / BITS_PER_WORD;
        return m_data[wIdx];
    }

    /**
     * Returns a bitmask for the desired bit within the word containing \a idx
     * \param idx The index of the desired bit
     */
    Word mask(Index idx) const
    {
        const Index wOffset = idx % BITS_PER_WORD;
        return Word(1) << wOffset;
    }

    /**
     * Returns a bitmask for the bits in the last word of the bitset.
     *
     * This function is only valid when isLasteWordFull() is false
     * \sa isLastWordFull()
     */
    Word lastWordMask() const
    {
        return mask(m_numBits ) - 1;
    }

    /**
     * Checks if index \idx corresponds to a valid index
     *
     * \note This function is a no-op in Release builds
     */
    void checkValidIndex(Index idx) const
    {
        AXOM_DEBUG_VAR(idx);
        SLIC_ASSERT_MSG(
            idx >= 0 && idx < m_numBits,
            "slam::Bitset attempted to out of range bit "
            << idx << ". Valid range is [0, " << m_numBits << ")."
        );
    }

    /**
     * Predicate to determine if we need special processing
     * for the final word of the bitset
     *
     * The last word is full when the bitset has exactly m_words * BITS_PER_WORD bits
     */
    bool isLastWordFull() const 
    { 
        const int lg = (1 << (LG_BITS_PER_WORD))-1;
        return (m_numBits & lg) == 0;
    }

private:
    ArrayType m_data;

    int m_numBits;
    int m_numWords;
};


// BitSet implementation here

BitSet::Index const BitSet::npos = -2;

void BitSet::clear()
{
    for (int i = 0; i < m_numWords; ++i)
    {
        m_data[i] = 0;
    }
}

void BitSet::set()
{
    if (m_numBits == 0) { return; }

    const Word ones = ~Word(0);
    for (int i = 0; i < m_numWords - 1; ++i)
    {
        m_data[i] = ones;
    }

    // Handle last word
    m_data[m_numWords - 1]
        = isLastWordFull() 
            ? ones
            : lastWordMask();
}

void BitSet::flip()
{
    if (m_numBits == 0) { return; }

    const Word ones = ~Word(0);
    for (int i = 0; i < m_numWords - 1; ++i)
    {
        m_data[i] ^= ones;
    }

    // Handle last word
    m_data[m_numWords - 1] 
        ^= isLastWordFull()
            ? ones
            : lastWordMask();
}

int BitSet::count() const
{
    int ctr = 0;

    for (int i = 0; i < m_numWords; ++i)
    {
        ctr += internal::popCount(m_data[i]);
    }
    return ctr;
}

bool BitSet::isValid() const
{
    bool valid = true;

    if (m_numBits < 0 || m_numWords < 0)
        valid = false;

    if (m_numBits == 0)
    {
        if (m_numWords != 1 || m_data[0] != Word(0))
        {
            valid = false;
        }
    }
    else
    {
        // check num words vs. num bits
        int expWords = (m_numBits - 1) / BITS_PER_WORD + 1;
        if (expWords != m_numWords)
            valid = false;

        // check that highest bits are not set
        if (!isLastWordFull())
        {
            const Word& lastWord = getWord(m_numBits, false);

            Word upperMask = ~lastWordMask();
            if ((lastWord & upperMask) != 0)
            {
                valid = false;
            }
        }
    }

    return valid;
}

BitSet::Index BitSet::find_first() const
{
    return find_next(-1);
}

BitSet::Index BitSet::find_next(Index idx) const
{
    // Handle boundary cases
    if (idx == npos || (idx+1 >= m_numBits))
    {
        return npos;
    }

    Index startWordIdx = 0;
    if (idx >= 0)
    {
        checkValidIndex(idx);

        const Index startIdx = idx + 1;
        startWordIdx = startIdx / BITS_PER_WORD;

        // Check for next set bit in current word
        const Index startOffset = startIdx % BITS_PER_WORD;

        const Word startWord = m_data[startWordIdx] >> startOffset;
        if (startWord != Word(0))
        {
            return (startWordIdx * BITS_PER_WORD)
                + internal::trailingZeros(startWord << startOffset);
        }

        ++startWordIdx;
    }

    // If not in current word, check remaining words
    for (int i = startWordIdx; i < m_numWords; ++i)
    {
        const Word& w = m_data[i];
        if (w != Word(0))
        {
            return (i * BITS_PER_WORD) 
                + internal::trailingZeros(w);
        }
    }
    return BitSet::npos;
}

BitSet& BitSet::operator|=(const BitSet& other)
{
    SLIC_ASSERT_MSG(size() == other.size(),
        "slam::BitSet Sizes must be the same for bit set operators."
        << " In operator|=(), BitSet has size " << size()
        << " and other BitSet has size " << other.size() << ".");

    for (int i = 0; i < m_numWords; ++i)
    {
        m_data[i] |= other.m_data[i];
    }

    return *this;
}

BitSet& BitSet::operator&=(const BitSet& other)
{
    SLIC_ASSERT_MSG(size() == other.size(),
        "slam::BitSet Sizes must be the same for bit set operators."
        << " In operator&=(), BitSet has size " << size()
        << " and other BitSet has size " << other.size() << ".");

    for (int i = 0; i < m_numWords; ++i)
    {
        m_data[i] &= other.m_data[i];
    }

    return *this;
}

BitSet& BitSet::operator^=(const BitSet& other)
{
    SLIC_ASSERT_MSG(size() == other.size(),
        "slam::BitSet Sizes must be the same for bit set operators."
        << " In operator^=(), BitSet has size " << size()
        << " and other BitSet has size " << other.size() << ".");

    for (int i = 0; i < m_numWords; ++i)
    {
        m_data[i] ^= other.m_data[i];
    }

    return *this;
}

BitSet& BitSet::operator-=(const BitSet& other)
{
    SLIC_ASSERT_MSG(size() == other.size(),
        "slam::BitSet Sizes must be the same for bit set operators."
        << " In operator-=(), BitSet has size " << size()
        << " and other BitSet has size " << other.size() << ".");

    for (int i = 0; i < m_numWords; ++i)
    {
        m_data[i] &= ~other.m_data[i];
    }

    return *this;
}

BitSet operator|(const BitSet & lhs, const BitSet & rhs)
{
    BitSet s(lhs);
    s |= rhs;
    return s;
}

BitSet operator&(const BitSet & lhs, const BitSet & rhs)
{
    BitSet s(lhs);
    s &= rhs;
    return s;
}

BitSet operator^(const BitSet & lhs, const BitSet & rhs)
{
    BitSet s(lhs);
    s ^= rhs;
    return s;
}

BitSet operator-(const BitSet & lhs, const BitSet & rhs)
{
    BitSet s(lhs);
    s -= rhs;
    return s;
}

bool BitSet::operator==(const BitSet & other) const
{
    if (size() != other.size() || m_numBits != other.m_numBits)
    {
        return false;
    }

    return m_data == other.m_data;
}

} // end namespace slam
} // end namespace axom

#endif //  SLAM_BITSET_H_
