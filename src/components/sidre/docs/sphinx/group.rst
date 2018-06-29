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
Group
==========

The Group class provides a tree structure.  Each Group has a name, a parent Group
(except for the root), zero or more child Groups, and zero or more Views.

Through a parent Group, Views and child Groups can be created, retrieved by name
or path, or iterated over.  A program can also query a parent Group about how many
child Groups and Views it has or the presence of a Group or View by
name.  Unlike many filesystems, the path may not contain the "parent" entry
(such as ".." on Unix filesystems).

Methods on the Group class use
`Conduit <https://github.com/LLNL/conduit>`_ to
:ref:`write the data (sub)tree <sidre-serial-io>` rooted in a Group to a file,
`HDF5 <https://www.hdfgroup.org/HDF5/>`_ handle, or other
Conduit protocol, or to an in-memory Conduit data structure.  The program may
provide an Attribute to the method call, so only Views with that Attribute
explicitly set will be written.  Other methods on the Group class allow loading
hierarchical data, constructing a tree structure rooted in the Group.


