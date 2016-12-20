.. _communicator_class_label:

Communicator Class
==================

The Communicator class is an abstract base class that defines the interface for
all Communicator classes.  Concrete instance need to inherit from this class and
implement these functions to be used when the Lumberjack class does any communication
work.

Functions
#########

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
##################

BinaryTreeCommunicator
********************

.. note:: This is the recommended Communicator.

This Communicator uses a standard Binary Tree design to scalably pass Messages between nodes.
Rank 0 is the root of the Binary Tree and the only node allowed to output messages.

RootCommunicator
********************

.. note:: This is not a recommended Communicator for production. It is provided for its simplistic design for debugging purposes.

This Communicator has all nodes directly connecting to the root node which
is rank 0.  The root node is the only node allowed to output messages.
