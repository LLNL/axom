******************************************************
Core concepts
******************************************************

Describe SPIO concept

IOManager
------

The API for the SPIO component is provided by class IOManager.  An IOManager
is constructed with an MPI communicator and does I/O operations on all ranks
associated with that communicator

The core functionality of IOManager is contained in the write and read methods.

  void write(sidre::DataGroup * group,
             int num_files,
             const std::string& file_string,
             const std::string& protocol);

write() is called in parallel with each rank passing in a DataGroup pointer
for its local data.  The calling code specifies the number of output files,
and IOManager organizes the output so that each file receives data from a
roughly equal number of ranks.  The files containing the data from the group
will have names of the format "file_string/file_string_*******.suffix", with a
7-digin integer value identifying the files from 0 to num_files-1, and the
suffix indicating the file format according tot he protocol argument.
Additionally write() will produce a root file with the name file_string.root
that holds some bookkeeping data about the other files and can also receive
extra user-specified data.

  void read(sidre::DataGroup * group,
            const std::string& root_file);

read() is called in parallel with the root file as an argument.  It must be
called on a run with the same processor count as the run that called write().
The first argument is a pointer to a group that contains no child groups or
vies, and the information in the root file is used to identify the files that
each processor will read to load data into the argument group.

The write and read methods above are sufficient to do a restart save/load
when the data is the group is completely owned by the Sidre data structures.
If Sidre is used to manage data that is externally allocated, the loading
procedure requires some additional restore data in the same externally-
allocated state.

First the read method is called, and the full hierarchy structure of the
group is loaded into the DataGroup, but no data is allocated for views
identified as external.  Then the calling code can examine the group and
allocate data for the external views.  DataView::setExternalDataPtr is used
to associate the pointer with the view.  Once this is done, IOManager's
loadExternalData can be used to load the data from the file into the
user-allocated arrays.

  code example here


Root file

  User can write additional data to root file
  blueprint index to provide metadata to VisIt

