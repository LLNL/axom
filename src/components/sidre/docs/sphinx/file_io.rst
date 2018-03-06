******************************************************
File I/O
******************************************************

Sidre provides for file I/O to HDF5 files and file handles, and to JSON files.
Serial I/O is accomplished using methods of the Group class; parallel I/O is done
using the IOManager class.

HDF5 is an optional dependency for Sidre.  Sidre APIs that rely on HDF5 will be
available only if Sidre was compiled with HDF5 support.

Serial File I/O
---------------

The Group class contains the :code:`save()`, :code:`load()`, and 
:code:`loadExternalData()` methods.  Each method can be called with a file 
name or an HDF5 handle.  The :code:`save()` and :code:`load()` methods 
allow the code to specify a protocol, specifying how the
operation should be performed and what file format should be used.  The
:code:`loadExternalData()` method takes no protocol argument since it is only used with
the sidre_hdf5 protocol.  The :code:`save()` method can also take a pointer to an
Attribute object.  If that Attribute pointer is null, all Views are saved, but
if an Attribute pointer is provided, only Views with that Attribute explicitly
set will be saved.

Insert a link to Doxygen for Group, and mention that protocols are listed here.

The :code:`load()` method retrieves the hierarchical data structure stored in the file 
or pointer, and instantiates Groups, Views and Attributes to represent the 
structure rooted in this Group.  By default, the contents of this Group are 
destroyed prior to reading the file's contents.  This can be suppressed by 
passing true as the preserve_contents argument to :code:`load()`, resulting in the 
current Group's subtree being merged with the file's contents.

Usage of :code:`save()` and :code:`load()` is shown in the following example.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: _serial_io_save_start
   :end-before: _serial_io_save_end
   :language: C++

The :code:`loadExternalData()` method is used to read "external" data from an HDF5 file
created with the sidre_hdf5 protocol.  This is data referred to by Views that is
not stored in Sidre Buffers but in a raw pointer.

The overloads of :code:`save()` and :code:`load()` that take HDF5 handles and the
:code:`loadExternalData()` method are public APIs, but also exist to implement parallel
I/O using the IOManager class.

Parallel File I/O
-----------------

To accomplish parallel I/O, Sidre provides the IOManager class.  This class is
instantiated with an MPI communicator and provides several overloads of the
:code:`write()` and :code:`read()` methods.  The :code:`IOManager::write()` 
method takes a Group
pointer, the number of files to write (distributed over MPI ranks), a string
specifying the base file name, and an I/O protocol (corresponding to the
protocols available to :code:`Group::save()` and :code:`load()`).  The 
:code:`IOManager::read()` method
takes a Group pointer into which it should read the data, a base file name, an
I/O protocol, and a boolean flag specifying if the contents of the Group should
be preserved.  By default, the Group is cleared before the data set is read in.

When :code:`write()` is called, the IOManager produces a root file and a number of data
files.  These files together contain the data hierarchy rooted in the Group
passed to :code:`write()`.  The :code:`read()` method will access all these files 
to instantiate the hierarchy under the Group pointer passed to it.

The IOManager constructor takes an optional boolean argument to use the SCR
library for scalable I/O management (i.e. use burst buffers if available).  This
defaults to false, so Sidre will use burst buffers only if specifically requested.  
IOManager provides other methods to add additional Groups or
Views to a root file, specified by name.

In the following example, an IOManager is created and used to write the contents
of the Group "root" in parallel.

.. literalinclude:: ../../tests/spio/spio_parallel.cpp
   :start-after: _parallel_io_save_start
   :end-before: _parallel_io_save_end
   :language: C++

