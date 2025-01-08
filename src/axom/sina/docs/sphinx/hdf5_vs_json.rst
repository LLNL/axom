.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _curvesets-label:

==========
HDF5 and JSON
==========

Sina's C++ and Python Code supports the formation of both HDF5 and JSON files
when exporting document objects.

Both options offer differing strengths and weaknesses when compared to each other,
so we've compiled data on the size, speed, and append efficiency of each file type
to give you a better idea of which file type is right for you.

==========
Why You Should Use JSON
==========

* JSON is more universally parsable
* JSON code is human readable allowing you more, easier flexibility with searching
for relevent data.  This is most applicable with smaller files since large amounts
of data or curve sets quickly become more efficient to navigate digitally
* JSON outperforms HDF5 in speed and size efficiency when dealing with non-curve set data 
and outperforms at smaller curve set sizes (just before 10^2 Curve Sets for Size and
around 10^3.25 for Speed)
* JSON offers more consistent append speeds, however on average these speeds are slower

==========
Why You Should Use HDF5
==========

* HDF5 offers better size and speed efficiency when dealing with larger files/curve sets 
and only outperforms more dramatically as size increases.
* Hierarchal structure leads to being 2.5x more size efficient and 5x faster at our largest
tested files
* HDF5 on the whole offers faster append times.  While its performance is far more subject to
variance than JSON its average speeds outperform JSON at all tested sizes.

==========
Images
==========

JSON vs HDF5 Size
--------


JSON vs HDF5 Speed
--------


JSON vs HDF5 Append
--------