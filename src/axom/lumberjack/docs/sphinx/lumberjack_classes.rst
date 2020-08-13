
Lumberjack Classes
==================

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

.. toctree::
   :maxdepth: 1

   lumberjack_class
   message_class
   communicator_class
   combiner_class
