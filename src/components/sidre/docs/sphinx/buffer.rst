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

.. _buffer-label:

==========
Buffer
==========

The Buffer class holds a linear data array of a specified length and type.
Buffers are allocated and managed through the containing Datastore class.  The
data stored in a Buffer is accessed through a View object or through the Buffer directly.


