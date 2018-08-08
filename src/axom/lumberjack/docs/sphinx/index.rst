
Lumberjack User Documentation
=============================

Lumberjack, named because it cuts down logs, is a C++ library that
provides scalable logging while reducing the amount of messages
written out the screen or file system.


Introduction
------------

Lumberjack was created to provide scalable logging with a simple programming
model while allowing developers to customize its behavior. It is named Lumberjack
because it cuts down logs. It uses MPI and a scalable binary tree reduction
scheme to combine duplicate messages and limit output to only the root node.


Requirements
------------

* MPI - MPI is fundamental to Lumberjack and without MPI is not useful.
* (Optional) C++11 - Can be optionally compiled with C++11 but not required.


Code Guarding
-------------

You tell if the ASC Toolkit was built with Lumberjack by the following
include and compiler define:

.. code-block:: c

    #include "axom/config.hpp"
    #ifdef ATK_USE_LUMBERJACK
        // Lumberjack work
    #endif


Classes
-------

Basic
*****

* :ref:`Lumberjack <lumberjack_class_label>` - Performs all high level functionality for the Lumberjack library.
* :ref:`Message <message_class_label>` - Holds all information pertaining to a Message.


Communicators
*************

Handles all node-to-node Message passing.

* :ref:`Communicator <communicator_class_label>` - Abstract base class that all Communicators must inherit from.
* :ref:`BinaryTreeCommunicator <binarytreecommunicator_class_label>` - Main Communicator that is implemented with a scalable Binary Tree scheme
* :ref:`RootCommunicator <rootcommunicator_class_label>` - non-scalable communication scheme that all nodes connect to the root node.  This is given for diagnostic purposes only.


Combiners
*********

Handles Message combination and tests whether Message classes should be combined.

* :ref:`Combiner <combiner_class_label>` - Abstract base class that all Combiners must inherit from.
* :ref:`TextEqualityCombiner <textequalitycombiner_class_label>` - Combines Message classes that have equal Text member variables.

**Contents:**

.. toctree::
   :maxdepth: 1

   core_concepts
   quick_start
   lumberjack_class
   message_class
   communicator_class
   combiner_class

