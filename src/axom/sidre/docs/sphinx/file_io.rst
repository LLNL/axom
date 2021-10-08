.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

******************************************************
Serial File I/O
******************************************************

Sidre provides file I/O to HDF5 files and file handles and to JSON files.
Serial I/O is accomplished using methods of the ``Group`` class; parallel I/O 
is done using the ``IOManager`` class.

HDF5 is an optional dependency for Sidre. Sidre APIs that rely on HDF5 will be
available only if Sidre was compiled with HDF5 support. The CMake variable
that controls dependency on HDF5 is ``AXOM_USE_HDF5``, defined in the
file axom/config.hpp.

.. _sidre-serial-io:

File I/O using the Group class
-------------------------------

The ``Group`` class provides methods: ``save()``, ``load()``, and
``loadExternalData()`` methods. Each method can be called with a file
name or an HDF5 handle.  The ``save()`` and ``load()`` methods allow users 
to specify a protocol, which indicates how the operation should be performed 
and what file format to use. The ``loadExternalData()`` method takes no 
protocol argument since it is only used with the ``sidre_hdf5`` protocol. 
The ``save()`` method can also take a pointer to an ``Attribute`` object. 
If the attribute pointer is null, all views are saved. If a non-null attribute 
pointer is given, only views with that attribute set will be saved.

The ``load()`` method retrieves the hierarchical data structure stored in the 
file and creates groups, views, and attributes to represent the hierarchy
with its root in the group on which the method is called. By default, the 
contents of the group are destroyed prior to reading the file contents. This 
can be suppressed by passing a Boolean value of true as the 
``preserve_contents`` argument to the ``load()`` method, resulting in the
current group subtree being merged with the file contents.

Usage of the ``save()`` and ``load()`` methods is shown in the following example.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _serial_io_save_start
   :end-before: _serial_io_save_end
   :language: C++

The ``loadExternalData()`` method is used to read "external" data from an HDF5 
file created with the ``sidre_hdf5`` protocol.  This is data referred to by 
*external* views. Such views refer to data that is not stored in Sidre buffers 
but in an external data allocation.

Overloads of the ``save()`` and ``load()`` methods that take HDF5 handles and 
the ``loadExternalData()`` method are used to implement parallel I/O through 
the Sidre ``IOManager`` class.

Please see `Sidre API Documentation <../../../../doxygen/html/sidretop.html>`_ 
for more detailed description of these ``Group`` class methods and the Sidre 
I/O protocols that can be used with them.

