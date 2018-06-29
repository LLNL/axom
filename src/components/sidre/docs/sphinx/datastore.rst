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

==========
Datastore
==========

The Datastore class provides a container for all Sidre data.  Each Datastore
contains one root Group, which provides access to the hierarchy of Group and View objects in
the Datastore.  A Datastore supports creation and destruction of Buffers and access
to individual Buffers via buffer IDs.
A Datastore also allows management of a list of Attributes
that can be set on any of its Views.


