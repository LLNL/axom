.. _core_concepts_label:

Core Concepts
=============

The following are core concepts required to understand how Lumberjack works.


Combining
---------

Combining Messages is how Lumberjack cuts down on the number of Messages output from
your program.  It does so by the giving the currently held Messages at the current node,
two at a time, to the Combiner classes that are currently registered to Lumberjack
when a Push happens.

Lumberjack only provides one Combiner, the TextEqualityCombiner. You can write your own
Combiners and register them with Lumberjack.  The idea is that each Combiner would have
its own criteria for whether a Message should be combined and how to combine that specific
Message with another of the same type.

Combiner's have two main functions, shouldMessagesBeCombined and combine.

The function shouldMessagesBeCombined, returns True if the criteria for has been
met that the given Messages should be combined.  For example in the TextEqualityCombiner,
if the Text strings are exactly equal, it signals they should be combined.

The function combine, takes two Messages and combines them in the way that is specific
to that Combiner class.  For example in the TextEqualityCombiner, the only thing
that happens is the second Message's ranks gets added to the first.  This is because
the text strings were equal.  This may not be the case for all Combiners that you write
yourself.


Pushing
-------

A push has three steps: combining, sending, and recieving Messages. When you queue
a Message into Lumberjack, it is held at the node that generated the Message until
you indicate the Lumberjack to push, either once or fully.  If you do not push,
then only the Messages generated at the root node will be outputed.

In a single push, nodes send their currently held Messages to the nodes their are connected
to based on the Communicator's communcation scheme.  For example in the BinaryTreeCommunicator,
children nodes send their Messages to their parent. While the root node only recieves Messages.
After a single push, it is not guaranteed that all Messages will be ready to be outputted.

A full push is a number of single pushes until all currently held Messages.  The Communicator
tells the Lumberjack class how many single pushes it takes to fully flush the system of
Messages.  For example in the BinaryTreeCommunicator, it is the log of the number of nodes.
