******************************************************
File I/O
******************************************************

Sidre provides for file I/O to HDF5 files and file handles and to JSON files.
Serial I/O is accomplished using methods of the Group class; parallel I/O is done
through the IOManager class.

HDF5 is an optional dependency for Sidre.  Sidre APIs that rely on HDF5 will be
available only if Sidre was compiled with HDF5 support.

Serial File I/O
---------------

The Group class contains the save(), load(), and loadExternalData() methods.
Each method can be called with a file name or an HDF5 handle.  The save() and
load() methods allow the user to specify a protocol, specifying how the
operation should be performed and what file format should be used.  The
loadExternalData() method takes no protocol argument since it is only used with
the sidre_hdf5 protocol.  The save() method can also take a pointer to an
Attribute object.  If that Attribute pointer is null, all Views are saved, but
if an Attribute pointer is provided, only Views with that Attribute explicitly
set will be saved.

Insert a link to Doxygen for Group, and mention that protocols are listed here.

The load() method retrieves the hierarchical data structure stored in the file 
or pointer, and instantiates Groups, Views and Attributes to represent the 
structure rooted in this Group.  By default, the contents of this Group are 
destroyed prior to reading the file's contents.  This can be suppressed by 
passing true as the preserve_contents argument to load(), resulting in the 
current Group's subtree being merged with the file's contents.

Usage of save() and load() is shown in the following example.

.. literalinclude:: ../../examples/sidre_createdatastore.cpp
   :start-after: serial_io_save_start
   :end-before: serial_io_save_end
   :language: C++

The loadExternalData() method is used to read "external" data from an HDF5 file
created with the sidre_hdf5 protocol.  This is data referred to by Views that is
not stored in Sidre Buffers but in a raw pointer.

The overloads of save() and load() that take HDF5 handles and the
loadExternalData() method are public APIs, but also exist to implement parallel
I/O using the IOManager class.

Parallel File I/O
-----------------

To accomplish parallel I/O, Sidre provides the IOManager class.  This class is
instantiated with an MPI communicator and provides several overloads of the
write() and read() methods.  The IOManager::write() method takes a Group
pointer, a number of files to write (distributed over MPI ranks), a string
specifying the base file name, and an I/O protocol (corresponding to the
protocols available to Group::save() and load()).  The IOManager::read() method
takes a Group pointer into which it should read the data, a base file name, an
I/O protocol, and a boolean flag specifying if the contents of the Group should
be preserved.  By default, the Group is cleared before the data set is read in.

When write() is called, the IOManager produces a root file and a number of data
files.  These files together contain the data hierarchy rooted in the Group
passed to write().  The read() method will access all these files to instantiate
the hierarchy under the Group pointer passed to it.

The IOManager constructor takes an optional boolean argument to use the SCR
library for scalable I/O management (i.e. use burst buffers if available).  This
defaults to false.  IOManager provides other methods to add additional Groups or
Views to a root file, specified by name.

In the following example, an IOManager is created and used to write the contents
of the Group "root" in parallel.

.. literalinclude:: ../../tests/spio/spio_parallel.cpp
   :start-after: parallel_io_save_start
   :end-before: parallel_io_save_end
   :language: C++

