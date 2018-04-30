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
* \brief Represents a dynamic bitset
*
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

/*! Union operator for two bit sets */
BitSet operator|(const BitSet & lhs, const BitSet & rhs);

/*! Intersection operator for two bit sets */
BitSet operator&(const BitSet & lhs, const BitSet & rhs);

/*! Xor operator for two bit sets */
BitSet operator^(const BitSet & lhs, const BitSet & rhs);

/*! Set difference operator for two bit sets */
BitSet operator-(const BitSet & lhs, const BitSet & rhs);

/**
* \class BitSet
* \brief A bitset class
*
* The initial design requires the size of the bitset to be
* given at construction time. This can be fixed through
* policies in the future.
*
* We do not yet support initializing the bits in a bitset using
* a constructor.
*
* Question: Do we need an internal 'reference' class?
*/
class BitSet
{
public:
    typedef int Index;
    typedef axom::common::uint64 Word;

    // Use vector at first -- to be updated with policy
    typedef std::vector<Word> ArrayType;

    static const Index npos;

private:

    enum
    {
        BITS_PER_WORD = sizeof(Word) << 3,
        LG_BITS_PER_WORD = 6
    };

public:

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

    BitSet(const BitSet& other)
        : m_data(other.m_data),
          m_numBits(other.m_numBits),
          m_numWords(other.m_numWords)
    {}

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

    int size() const { return m_numBits; }

    // Return the number of set bits
    int count() const;
    
    /** BitSet union-asignment operator */
    BitSet& operator|=(const BitSet& other);

    /** BitSet interection-asignment operator */
    BitSet& operator&=(const BitSet& other);

    /** BitSet xor-asignment operator */
    BitSet& operator^=(const BitSet& other);

    /** BitSet difference-asignment operator */
    BitSet& operator-=(const BitSet& other);

    /*! Equality operator for two bit sets */
    bool operator==(const BitSet & other) const;

    /*! Inequality operator for two bit sets */
    bool operator!=(const BitSet & other) const
    {
        return !(operator==(other));
    }

    Index find_first() const;
    Index find_next(Index idx) const;
public:
    // Operations that affect all bits

    void clear();
    void set();
    void flip();

public:
    // Operations that affect a single bit

    void clear(Index idx) { getWord(idx) &= ~mask(idx); }
    void set(Index idx)   { getWord(idx) |= mask(idx);  }
    void flip(Index idx)  { getWord(idx) ^= mask(idx);  }

    bool test(Index idx)  const
    {
        return (getWord(idx) & mask(idx)) != Word(0);
    }

    bool isValid() const;

private:
    Word& getWord(Index idx, bool checkIndexValid = true)
    {
        if(checkIndexValid)
            checkValidIndex(idx);

        const Index wIdx = idx / BITS_PER_WORD;
        return m_data[wIdx];
    }

    const Word& getWord(Index idx, bool checkIndexValid = true) const
    {
        if (checkIndexValid)
            checkValidIndex(idx);

        const Index wIdx = idx / BITS_PER_WORD;
        return m_data[wIdx];
    }

    Word mask(Index idx) const
    {
        const Index wOffset = idx % BITS_PER_WORD;
        return Word(1) << wOffset;
    }

    Word lastWordMask() const
    {
        return mask(m_numBits ) - 1;
    }

    void checkValidIndex(Index idx) const
    {
        AXOM_DEBUG_VAR(idx);
        SLIC_ASSERT_MSG(
            idx >= 0 && idx < m_numBits,
            "slam::Bitset attempted to out of range bit "
            << idx << ". Valid range is [0, " << m_numBits << ")."
        );
    }

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
