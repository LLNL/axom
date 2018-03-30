.. _spio-core-concepts:
******************************************************
Core concepts
******************************************************

IOManager
---------

The API for the SPIO component is provided by class IOManager.  An IOManager
is constructed with an MPI communicator and does I/O operations on all ranks
associated with that communicator

The core functionality of IOManager is contained in the write and read methods.

.. code-block:: cpp

  void write(sidre::DataGroup * group,
             int num_files,
             const std::string& file_string,
             const std::string& protocol);

write() is called in parallel with each rank passing in a DataGroup pointer
for its local data.  The calling code specifies the number of output files,
and IOManager organizes the output so that each file receives data from a
roughly equal number of ranks.  The files containing the data from the group
will have names of the format "file_string/file_string_*******.suffix", with a
7-digit integer value identifying the files from 0 to num_files-1, and the
suffix indicating the file format according to the protocol argument.
Additionally write() will produce a root file with the name file_string.root
that holds some bookkeeping data about the other files and can also receive
extra user-specified data.

.. code-block:: cpp

  void read(sidre::DataGroup * group,
            const std::string& root_file);

read() is called in parallel with the root file as an argument.  It must be
called on a run with the same processor count as the run that called write().
The first argument is a pointer to a group that contains no child groups or
views, and the information in the root file is used to identify the files that
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

Below is a code example for loading external data.  We assume that this code
somehow has knowledge that the root group contains a single external view at
the location "fields/external_array" describing an array of doubles.  See the
Sidre documentation for information about how to query the Sidre data
structures for this type of information when the code does not have a priori
knowledge.

.. code-block:: cpp

  // Construct a DataStore with an empty root group.
  DataStore * ds = new DataStore();
  DataGroup * root = ds->getRoot();

  // Read from file into the root group.  The full Sidre hierarchy is built,
  // but the external view is created without allocating a data buffer.
  IOManager reader(MPI_COMM_WORLD);
  reader.read(root, "checkpoint.root");

  // Get a pointer to the external view. 
  DataView * external_view = root->getView("fields/external_array");

  // Allocate storage for the array and associate it with the view.
  double * external_array = new double[external_view->getNumElements()];
  external_view->setExternalDataPtr(external_array);

  // Load the data values from file into the external view.
  reader.loadExternalData(root, "checkpoint.root");


User-specified data in the root file
------------------------------------

The root file is automatically created to provide the IOManager with
bookkeeping information that is used when reading data, but it can also
be used to store additional data that may be useful to the calling code or
is needed to allow other tools to interact with the data in the output files,
such as for visualization.  For example, Conduit's blueprint index can be
stored in a DataGroup written to the root file to provide metadata about the
mesh layout and data fields that can be visualized from the output files.

.. code-block:: cpp

  void writeGroupToRootFile(sidre::DataGroup * group,
                            const std::string& file_name);

  void writeGroupToRootFileAtPath(sidre::DataGroup * group,
                                  const std::string& file_name,
                                  const std::string& group_path);

  void writeViewToRootFileAtPath(sidre::DataView * view,
                                 const std::string& file_name,
                                 const std::string& group_path);

The above methods are used to write this extra data to the root file.  The
first simply writes data from the given group to the top of the root file,
while the latter two methods write their Sidre objects to a path that must
already exist in the root file.
