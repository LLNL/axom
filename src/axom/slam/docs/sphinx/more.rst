.. ## Copyright (c) Lawrence Livermore National Security, LLC and other
.. ## Axom Project Contributors. See top-level LICENSE and COPYRIGHT
.. ## files for dates and other details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

====================
Additional utilities
====================

BitSet
======

``BitSet`` stores a runtime-sized sequence of bits in a packed array. Use it
to record a selection of mesh entities when you need one flag per position.
It supports union, intersection, difference and exclusive-or operations on
bitsets of equal size.

``BitSet`` does not provide ``operator[]``, so use ``test(i)`` to inspect a bit. 
To visit only the set bits, use ``find_first()`` and ``find_next(i)`` until
they return ``BitSet::npos``. ``count()`` reports the number of set bits,
while ``size()`` reports the total number of bit positions.

ModularInt
==========

``ModularInt`` wraps an integer into the range ``[0, modulus())``. It is useful
for circular traversal, such as moving to the next vertex of a polygon without
a separate test at the last vertex. Increment, decrement and arithmetic
normalize the stored value to that range.

Supply a positive modulus at construction, or fix it at compile time through
the size policy. Include ``axom/slam/ModularInt.hpp`` to use this utility.
