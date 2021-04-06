.. ## Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _group-label:

==========
Group
==========

Sidre Group objects are used to define a tree-like hierarchical organization
for application data, such as meshes and fields used in a simulation. Each 
Group has a name and one parent Group (except for the root Group) and contains
zero or more child Groups and zero or more Views. A Sidre DataStore has 
exactly one root Group which is created when the DataStore object
is created. The root Group's name is initially the empty string.
See :ref:`datastore-label` for more information.

A Group hierarchy is constructed by creating Groups that are children of the
root Group, children of those Groups, and so on. All Groups in a subtree 
rooted at a particular Group are considered descendants of that Group. View
objects can be created in Groups to hold or provide access to data.

.. note:: Group and View objects can only be created and destroyed using
          Group methods provided for this. The Group and View constructors 
          and destructors are private. 

A Group or View object is owned by the Group that created it; i.e., its parent 
Group or *owning* Group, respectively. Groups and Views maintain pointers to 
their parent/owning Group. Thus, one may *walk* up or down a Group hierarchy
to access different Groups and Views in it.

.. note:: * The name (string) of a Group or View **must be unique** within its 
            parent/owning Group.
          * A Group or View has a unique integer identifier within its
            parent/owning group, which is generated when it is created.
          * Views and child Groups in a Group can be accessed by name or
            integer id.

A Group object can be moved or copied to another Group. When a Group is moved
to another Group, it is removed from its original parent and the Group to 
which it is moved becomes its parent. This implies that the entire subtree 
of Groups and Views within the moved Group is moved as well and can no longer 
be accessed via the original parent Group. When a Group is copied to another 
Group, a copy of the entire Group subtree rooted at the copied Group is added
to the Group to which it is copied. A **shallow** copy is performed for the
data in each View; i.e., a new View object is created in the destination, but 
the data is shared by the original and new View.

.. note:: View copy operations perform **shallow** copies of the View data.

Some methods for creating, destroying, querying, and retrieving Groups and
Views take a string with *path syntax*, where parent and child Group names
are joined with the path separator character, '/'.
Other methods take the name of an
immediate child of a Group. Methods that require the name of a direct child 
are marked with 'Child' in their name, such as ``hasChildView()`` and 
``hasChildGroup()``. When a path string is passed to a method that accepts 
path syntax, the last item in the path indicates the item to be created, 
destroyed, accessed, etc.  For example,::

   View* view = group->createView("foo/bar/baz");

is equivalent to::

   View* view = group->createGroup("foo")->createGroup("bar")->createView("baz");

In particular, intermediate Groups "foo" and "bar" will be created in this 
case if they don't already exist. The path syntax is similar to a Unix 
filesystem, but the path string **may not** contain the parent entry
(such as "../foo").

----------------------------
Methods to Operate on Groups
----------------------------

The following lists summarize Group methods that support operations related to 
Group objects.

.. note:: * Methods that access Groups by index only work with the direct 
            children of the current Group because an id has no meaning 
            outside of the indexing of the current group. None of these methods 
            is marked with 'Child' in its name.
          * When Groups are created, destroyed, copied, or moved,
            ids of other Views and Groups in parent Group objects may
            become invalid. This is analogous to iterator invalidation for
            containers when the container contents change.

Create, Modify, and Destroy Groups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

 * Create a child Group given a name (child) or path (other descendant). 
   If a path is given, intermediate Groups in path are created, if needed. 
 * Rename a Group.  A Group cannot be renamed to the empty string, to
   a string containing the path separator character, or to the name of
   another Group or View owned by the same parent.
 * Destroy a descendant Group with given id (child), or name/path (child or 
   other descendant).
 * Destroy all child groups in a Group.

.. note:: When a Group is destroyed, all Groups and Views in the subtree 
          rooted at the destroyed Group are also destroyed. However, the 
          data associated with those Views will remain intact.

Group Properties  
^^^^^^^^^^^^^^^^^^^^^^^

 * Retrieve the name or id of a Group object.
 * Retrieve the full path name from the root of the tree to a Group object.
 * Get a pointer to the parent Group of a Group.
 * Query the number of child Groups of a Group.
 * Query whether a Group has a descendant Group with a given name or path.
 * Query whether a Group has a child Group with a given integer id.
 * Query the name of a child Group with a given id, or the id of a child Group
   with a given name.
 * Get a pointer to the DataStore that owns the hierarchy in which a Group 
   resides.

Group Access
^^^^^^^^^^^^

 * Retrieve an immediate child Group with a given name or id, or a descendant
   Group with a given path.
 * Iterate over the set of child Groups in a Group.

Move and Copy Groups
^^^^^^^^^^^^^^^^^^^^^^

 * Move a Group, and its associated subtree, from its parent Group and make it
   a child of another Group.
 * Create a copy of Group subtree rooted at some Group and make it a child of 
   another Group.
 * Query whether Group subtree is equivalent to another; i.e., identical 
   subtree structures with same names for all Groups and Views, and Views are 
   also equivalent (see :ref:`view-interface-label`).

----------------------------
Methods to Operate on Views
----------------------------

The Group methods that support operations related to View objects are 
summarized below. For more details on View concepts and operations, please
see :ref:`view-label`.

.. note:: Methods that access Views by index work only with the
          Views owned by the current Group because an id has no meaning 
          outside of the indexing of the current group. None of these methods 
          is marked with 'Child' in its name.

Create Views
^^^^^^^^^^^^^

 * Create a View in the Group with a name only.
 * Create a View in the Group with a name and data description.
 * Create a View in the Group with a name and with a Buffer attached. The
   View may or may not have a data description.
 * Create a View in the Group with a name and an external data pointer. The
   data may or may not be described.
 * Create a View in the Group with a name and data description, and allocate
   the data. Implicitly the data is held in a Buffer that is attached to the
   View.
 * Create a View in the Group with a name holding a given scalar or string.

Destroy Views
^^^^^^^^^^^^^^

 * Destroy View with given id (child), or name/path (View in the Group or some 
   descendant Group), and leave View data intact.
 * Destroy all Views in the Group, and leave their data intact.
 * Destroy View with given id, or name/path, and destroy their data.
 * Destroy all Views in the Group and destroy their data.

View Queries
^^^^^^^^^^^^^^^^

 * Query the number of Views in a Group.
 * Query whether a Group subtree has a View with a given name or path.
 * Query whether a Group has a View with a given integer id.
 * Query the name of a View with a given id, or the id of a View with a given 
   name.

View Access
^^^^^^^^^^^^^

 * Retrieve a View with a given name or id, or a descendant View (somewhere
   in the subtree) with a given path.
 * Iterate over the set of Views in a Group.

Move and Copy Views
^^^^^^^^^^^^^^^^^^^^

 * Move a View from its owning Group to another Group (removed from original 
   owning Group).
 * Copy a View to another Group. Note that this is **shallow** copy of the
   View data; i.e., it is shared by the original and new View.

----------------------------
Group I/O Operations
----------------------------

The Group interface provides methods to perform data I/O operations on Views
in the Group subtree rooted at any Group.

 * Copy a description of a Group subtree to a conduit::Node.
 * Create native and external data layouts in conduit::Node hierarchies 
   (used mainly for I/O operations)
 * Save and load Group subtrees, including data in associated Views, to and
   from files. A variety of methods are provided to support different I/O
   operations, different I/O protocols, etc.

I/O methods on the Group class use `Conduit <https://github.com/LLNL/conduit>`_
to :ref:`write the data (sub)tree <sidre-serial-io>` rooted in a Group to a 
file, `HDF5 <https://www.hdfgroup.org/HDF5/>`_ handle, or other
Conduit protocol, or to an in-memory Conduit data structure. An application may
provide an Attribute to the method call, so only Views with that Attribute
explicitly set will be written. See :ref:`spio-core-concepts` for more 
information.

