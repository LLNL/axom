
SPIO User Documentation
=========================

The SPIO (Sidre Parallel Input/Ouptput) component of the CS Toolkit provides
an interface to manage the parallel I/O of the data managed by Sidre.  It
enables the writing of data from parallel runs and can
be used for the purposes of restart or visualization.

Introduction
-------------

SPIO relies on the fact that Sidre's DataGroup and DataView objects are
capable of saving and loading themselves.  These I/O operations in Sidre are
inherently serial, so SPIO provides the coordination of the I/O of multiple
Sidre objects that exists across all the MPI ranks of a parallel run.

* The internal details of the I/O of individual Sidre objects are opaque to
  SPIO, which needs only to make calls to the I/O methods in Sidre's public
  API.
* Sidre data is written from M ranks to N files (M >= N), and the files can
  be read to restart a run on M ranks.
* When saving output, a root file is created that contains some bookkeeping
  data that is used to coordinate a subsequent restart read.
* The calling code can also add extra data to the root file to provide 
  metadata that gives necessary instructions to visualization tools.


**Contents:**

.. toctree::
   :maxdepth: 2

   first_example
   core_concepts

`Class documentation <../../doxygen/html/annotated.html>`_

