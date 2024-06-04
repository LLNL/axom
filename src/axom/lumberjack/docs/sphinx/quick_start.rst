.. _quick_start_label:

Quick Start
===========

This quick start guide goes over the bare minimum you need to do to get up
and running with Lumberjack.  You can find this example in the repository under
Lumberjack's examples directory.

This example uses the Binary Tree Communicator and queues one unique message and
three similar messages per rank.  They are combined and then :ref:`pushed <push_label>`
fully through the tree.

The following files need to be included for Lumberjack:

.. code-block:: c

    # Lumberjack specific header
    #include "axom/lumberjack.hpp"

    # MPI and C++
    #include <mpi.h>
    #include <iostream>


Basic MPI setup and information:

.. code-block:: c

    // Initialize MPI and get rank and comm size
    MPI_Init(&argc, &argv);

    int commRank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
    int commSize = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &commSize);


Initialize Lumberjack:

.. code-block:: c

    // Determine how many ranks we want to individually track per message
    int ranksLimit = commSize/2;

    // Initialize which lumberjack communicator we want
    axom::lumberjack::BinaryTreeCommunicator communicator;
    communicator.initialize(MPI_COMM_WORLD, ranksLimit);

    // Initialize lumberjack
    axom::lumberjack::Lumberjack lj;
    lj.initialize(&communicator, ranksLimit);


This queues the individual messages into Lumberjack:

.. code-block:: c

    // Queue messages into lumberjack
    if (commRank == 0){
        lj.queueMessage("This message will not be combined");
    }
    else {
        lj.queueMessage("This message will be combined");
        lj.queueMessage("This message will be combined");
        lj.queueMessage("This message will be combined");
    }


This is how you fully :ref:`push <push_label>` all Messages through the Communicator,
which also :ref:`combines <combine_label>` Messages before and after pushing :

.. code-block:: c

    // Push messages fully through lumberjack's communicator
    lj.pushMessagesFully();


Optionally, you could spread the :ref:`pushing <push_label>` over the
course of your work by doing the following:

.. code-block:: c

    int cycleCount = 0;
    int cycleLimit = 10;
    for (int i = 0; i < someLoopLength; ++i){
        //
        // Do some work
        //
        lj.queueMessage("This message will combine")
        ++cycleCount;
        if (cycleCount > cycleLimit) {
            // Incrementally push messages through system
            lj.pushMessagesOnce();
            cycleCount = 0;
        }
    }


Once you are ready to retrieve your messages, do so by the following:

.. code-block:: c

    // Determine if this is an output node
    if (lj.isOutputNode()){
        // Get Messages from Lumberjack
        std::vector<axom::lumberjack::Message*> messages = lj.getMessages();
        for(int i=0; i<(int)(messages.size()); ++i){
            // Output a single Message at a time to screen
            std::cout << "(" << messages[i]->stringOfRanks() << ") " << messages[i]->count() <<
                         " '" << messages[i]->text() << "'" << std::endl;
        }
        // Clear already outputted Messages from Lumberjack
        lj.clearMessages();
    }

Finalize Lumberjack, the Lumberjack Communicator and MPI in the following order to guarantee nothing
goes wrong:

.. code-block:: c

    // Finalize lumberjack
    lj.finalize();
    // Finalize the lumberjack communicator
    communicator.finalize();
    // Finalize MPI
    MPI_Finalize();

