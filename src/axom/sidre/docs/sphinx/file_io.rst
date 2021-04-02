.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Serial File I/O
******************************************************

Sidre provides for file I/O to HDF5 files and file handles and to JSON files.
Serial I/O is accomplished using methods of the Group class; parallel I/O is done
using the IOManager class.

HDF5 is an optional dependency for Sidre.  Sidre APIs that rely on HDF5 will be
available only if Sidre was compiled with HDF5 support.  The symbol that controls
dependency on HDF5 is ``AXOM_USE_HDF5``, defined in the source file axom/config.hpp.

.. _sidre-serial-io:

File I/O using Group class
--------------------------

The Group class contains the ``save()``, ``load()``, and
``loadExternalData()`` methods.  Each method can be called with a file
name or an HDF5 handle.  The ``save()`` and ``load()`` methods
allow the code to specify a protocol, specifying how the
operation should be performed and what file format should be used.  The
``loadExternalData()`` method takes no protocol argument since it is only used with
the sidre_hdf5 protocol.  The ``save()`` method can also take a pointer to an
Attribute object.  If that Attribute pointer is null, all Views are saved, but
if an Attribute pointer is provided, only Views with that Attribute explicitly
set will be saved.

Insert a link to Doxygen for Group, and mention that protocols are listed here.

The ``load()`` method retrieves the hierarchical data structure stored in the file
or pointer, and instantiates Groups, Views and Attributes to represent the
structure rooted in this Group.  By default, the contents of this Group are
destroyed prior to reading the file's contents.  This can be suppressed by
passing true as the preserve_contents argument to ``load()``, resulting in the
current Group's subtree being merged with the file's contents.

Usage of ``save()`` and ``load()`` is shown in the following example.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _serial_io_save_start
   :end-before: _serial_io_save_end
   :language: C++

The ``loadExternalData()`` method is used to read "external" data from an HDF5 file
created with the sidre_hdf5 protocol.  This is data referred to by Views that is
not stored in Sidre Buffers but in a raw pointer.

The overloads of ``save()`` and ``load()`` that take HDF5 handles and the
``loadExternalData()`` method are public APIs that are used to implement parallel
I/O through the IOManager class.

