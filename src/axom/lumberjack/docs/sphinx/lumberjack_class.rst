.. _lumberjack_class_label:

Lumberjack Class
================

The Lumberjack class is where all high-level functionality of the library is done,
such as adding, retrieving, and combining messages and telling the given Communicator
to push Messages through the communication scheme.  You can also add and remove
Combiner classes, as well as tell if the current node is supposed to output any messages.


Functions
---------

General
^^^^^^^

============== ===================
Name           Description
============== ===================
initialize     Starts up Lumberjack. Must be called before anything else.
finalize       Cleans up Lumberjack. Must be called when done with Lumberjack.
isOutputNode   Returns whether this node should output messages.
ranksLimit     Sets the limit on individually tracked ranks
ranksLimit     Gets the limit on individually tracked ranks
============== ===================

Combiners
^^^^^^^^^

============== ===================
Name           Description
============== ===================
addCombiner    Adds a combiner to Lumberjack
removeCombiner Removes a specific combiner from Lumberjack
clearCombiners Removes all currently registered Combiners from Lumberjack
============== ===================

Messages
^^^^^^^^

================== ===================
Name               Description
================== ===================
clearMessages      Delete all Messages currently held by this node.
getMessages        Get all Messages currently held by this node.
queueMessage       Adds a Message to Lumberjack
pushMessagesOnce   Moves Messages up the communication scheme once
pushMessagesFully  Moves all Messages through the communication scheme to the output node.
================== ===================

