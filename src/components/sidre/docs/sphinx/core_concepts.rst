******************************************************
Core concepts
******************************************************

Sidre provides four main classes: Datastore, Buffer, Group, and View.  In
combination, these classes implement a data store with a tree structure:

* DataStore acts as an overall container,
* Buffer holds the data,
* Group provides the hierarchical tree structure and I/O operations,
* View gives access to the data in the Buffers.  

Sidre also provides the Attribute class, allowing a program to attach metadata 
to View objects.

Datastore
---------

The Datastore class provides a container for all Sidre data.  Each Datastore
contains a root Group, which provides access to the Group and View objects of
the Datastore.  Datastores allow allocation, querying, iteration over, and
disposal of Buffers.  A Datastore also allows management of a list of Attributes
that can be set on any of its Views.

Buffer
------

The Buffer class holds a linear data array of a specified length and type.
Buffers are allocated and managed through the containing Datastore class.  The
data stored in a Buffer is accessed through a View object.

Group
------

The Group class provides a tree structure.  Each Group has a name, a parent Group
(except for the root), zero or more child Groups, and zero or more child Views.

Through a parent Group, child Groups and Views can be created, retrieved by name
or path, or iterated over.  A program can also query a parent Group about how many
child Groups and Views it has or the presence of a Group or View by
name.  Unlike many filesystems, the path may not contain the "parent" entry
(such as ".." on Unix filesystems).

Methods on the Group class (:ref:`discussed here <sidre-serial-io>`) use
`Conduit <https://github.com/LLNL/conduit>`_ to
write the data (sub)tree rooted in a Group to a file,
`HDF5 <https://www.hdfgroup.org/HDF5/>`_ handle, or other
Conduit protocol, or to an in-memory Conduit data structure.  The program may
provide an Attribute to the method call, so only Views with that Attribute
explicitly set will be written.  Other methods on the Group class allow loading
hierarchical data, constructing a tree structure rooted in the Group.

View
------

The View class provides access to data.  A View has a parent Group, a name, a
data type, a length, an offset, and a stride; this information is available to
programs that use Sidre.  A program can obtain a pointer from the View and use
the pointer to retrieve and store data.  A View can refer to a Buffer or to an
external pointer.  If a code uses an external pointer, the View can store data 
type, length, offset, and stride supplied by the code for later use.  A View may 
also refer to an opaque data pointer, recording no information beyond the bare
pointer.

Attribute
---------

Attributes provide storage for View metadata.  When each
Attribute is constructed in the DataStore, it gets a name and a default value.
Each View inherits all of its DataStore's Attributes at their default values.
A program may explicitly set any of these Attributes for each
View.  The program may also query the value of a View's Attribute, query whether
the Attribute was explicitly set, or clear the Attribute back to its default
value.  

Attribute values are available for a program to use in its own logic.  If a
program provides an Attribute pointer to :code:`Group::save()` (discussed in the next
section), only Views with that Attribute explicitly set will be saved.  Further
extensions to Sidre that use Attributes and their values are planned.

.. Is this brief note about AttrValues needed?

The AttrValues class exists in Sidre as a helper to connect Attributes with
Views and provides no external functionality.

