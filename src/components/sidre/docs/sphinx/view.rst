.. ##
.. ## Copyright (c) 2017-18, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory
.. ##
.. ## LLNL-CODE-741217
.. ##
.. ## All rights reserved.
.. ##
.. ## This file is part of Axom.
.. ##
.. ## For details about use and distribution, please read axom/LICENSE.
.. ##

.. _view-label:

==========
View
==========

The View class provides applications with access to data pointers.  A View is
owned by one Group and has a name and a pointer to data.  A View interprets the
data with a data type, length (number of elements), offset, and stride; this
information is available to programs that use Sidre.  The View's pointer can
refer to data stored in a Buffer or to a memory location outside the Datastore.
For an external pointer, the View can store data type, length, offset, and
stride supplied by the code.  A View may also refer to an opaque data pointer,
recording no information beyond the bare pointer.

A View may also hold scalar data in the form of a string or number.
