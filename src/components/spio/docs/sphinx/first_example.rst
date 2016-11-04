******************************************************
An introductory example
******************************************************

Introductory example for SPIO.  First, an example of writing output.

.. code-block:: cpp

  // Get values for MPI rank and size
  int my_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  int num_ranks;
  MPI_Comm_size(MPI_COMM_WORLD, &num_ranks);

  // Determine number of output files desired.
  // num_ranks/2 is simple choice for example purposes.
  int num_output = num_ranks / 2;
  if (num_output == 0) {
    num_output = 1;
  }

  // Create a DataStore and get the root group
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  // Create a Group called "domain"
  Group * domain_group = root->createGroup("domain");

  // Code to fill domain_group with data omitted

  // Create IOManager and write 
  IOManager writer(MPI_COMM_WORLD);
  writer.write(domain_group, num_output, "domain_output", "sidre_hdf5");
 
Now, an example of reading files written by the above example

.. code-block:: cpp

  // Create a DataStore and get the root group
  DataStore * ds = new DataStore();
  Group * root = ds->getRoot();

  // Create a Group called "domain"
  Group * domain_group = root->createGroup("domain");

  // Create IOManager and read from files into domain_group.
  IOManager reader(MPI_COMM_WORLD);
  reader.read(domain_group, "domain_output.root");
 
