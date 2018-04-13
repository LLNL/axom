******************************************************
Core concepts
******************************************************

Sidre provides five main classes: Datastore, Buffer, Group, View, and Attribute.  In
combination, these classes implement a data store with a tree structure:

* DataStore acts as an overall container.  It is the main interface to access a data hierarchy.
* Buffer describes and holds data in memory.
* Group defines parent-child relationships in a hierarchical tree data structure provides serial I/O operations.
* View provides a description of data and gives access to it.
* Attribute allows a program to attach metadata to View objects.

Datastore
---------

The Datastore class provides a container for all Sidre data.  Each Datastore
contains one root Group, which provides access to the hierarchy of Group and View objects in
the Datastore.  A Datastore supports creation and destruction of Buffers and access
to individual Buffers via buffer IDs.
A Datastore also allows management of a list of Attributes
that can be set on any of its Views.

Buffer
------

The Buffer class holds a linear data array of a specified length and type.
Buffers are allocated and managed through the containing Datastore class.  The
data stored in a Buffer is accessed through a View object or through the Buffer directly.

Group
------

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

View
------

The View class provides applications with access to data pointers.  A View is
owned by one Group and has a name and a pointer to data.  A View interprets the
data with a data type, length (number of elements), offset, and stride; this
information is available to programs that use Sidre.  The View's pointer can
refer to data stored in a Buffer or to a memory location outside the Datastore.
For an external pointer, the View can store data type, length, offset, and
stride supplied by the code.  A View may also refer to an opaque data pointer,
recording no information beyond the bare pointer.

A View may also hold scalar data in the form of a string or number.

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

