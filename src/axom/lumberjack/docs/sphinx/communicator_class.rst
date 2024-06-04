.. _communicator_class_label:

Communicator Class
==================

The Communicator class is an abstract base class that defines the interface for
all Communicator classes.  Concrete instances need to inherit from this class and
implement these functions to be used when the Lumberjack class does any communication
work.

Functions
---------

========================= ===================
Name                      Description
========================= ===================
initialize                Starts up the Communicator. Must be called before anything else.
finalize                  Cleans up the Communicator. Must be called when finished.
rank                      Returns the rank of the current node.
ranksLimit                Getter/Setter for the limit on individually stored ranks.
numPushesToFlush          Returns the number of individual pushes to completely flush all Messages.
push                      Pushes all currently held Messages once up structure.
isOutputNode              Returns whether this node should output messages.
========================= ===================

Concrete Instances
------------------

.. _binarytreecommunicator_class_label:

BinaryTreeCommunicator
^^^^^^^^^^^^^^^^^^^^^^

.. note:: This is the recommended Communicator.

This Communicator uses a standard Binary Tree design to scalably pass Messages between nodes.
Rank 0 is the root of the Binary Tree and the only node allowed to output messages. For each single
push, the child nodes send their currently held messages to their parents without waiting to
receive messages themselves.  For a full push, this communicator takes the log of nodes to completely flush
all currently held messages to the root node.

.. _rootcommunicator_class_label:

RootCommunicator
^^^^^^^^^^^^^^^^

.. note:: This Communicator is useful for debugging purposes, but will not scale as well as the recommended BinaryTreeCommunicator.

This Communicator has all nodes directly connecting to the root node which
is rank 0.  The root node is the only node allowed to output messages.
Each single push, the child nodes send their currently held messages
to the root.  After each push the tree is completely flushed.
