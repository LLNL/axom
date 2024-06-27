.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _relationships-label:

=============
Relationships
=============

Sina ``Relationship`` objects are used as a way to describe the correlation between
two ``Record`` instances (and/or ``Record`` inheritors, e.g. ``Run``). A ``Relationship``
consists of three parts: a subject, an object, and a predicate. The subject and object
must be IDs referring to valid instances of ``Records``, while the predicate may be
any string.

In describing the connection between objects, a ``Relationship`` is read as "<subject>
<predicate> <object>". For example, in the relationship "Alice knows Bob", "Alice" is
the subject, "knows" is the predicate, and "Bob" is the object. Below are some additional
examples:

    - Task_22 contains Run_1024
    - msub_1_1 describes out_j_1_1
    - Carlos sends an email to Dani
    - local_task_12 runs before local_run_14

Note that a ``Relationship`` is described in the active voice. **Avoiding the passive
voice in predicates is recommended**, as this keeps the "direction" of the relationship
constant. An example of a passively-voiced ``Relationship`` is "Dani is emailed by Carlos".
Instead, this should be phrased as "Carlos emails Dani".

Below is an example showcasing how to construct a ``Relationship`` programmatically. Here,
we assemble a ``Relationship`` showing that "Task_22 contains Run_1024":

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create IDs for both Task 22 and Run 1024
        axom::sina::ID task22{"Task_22", sina::IDType::Global};
        axom::sina::ID run1024{"Run_1024", sina::IDType::Global};

        // Create the relationship and print it out
        axom::sina::Relationship myRelationship{task22, "contains", run1024};
        std::cout << myRelationship.toNode().to_json() << std::endl;
    }

If executed, the above code will output:

.. code:: json

    {
        "object": "Run_1024",
        "predicate": "contains",
        "subject": "Task_22"
    }

As with any other Sina ID, the subject or object may be either local (uniquely refer to
one object in a Sina file) or global (uniquely refer to one object in a database). Local
IDs are replaced with global ones upon ingestion; all ``Relationship`` instances referring
to that local ID (as well as the ``Record`` possessing that ID) will be updated to use the
same global ID.

Let's add on to our previous example to demonstrate this:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create IDs for Task 22 and Run 1024
        axom::sina::ID task22{"Task_22", sina::IDType::Global};
        axom::sina::ID run1024{"Run_1024", sina::IDType::Global};

        // Create the relationship and print it out
        axom::sina::Relationship myRelationship{task22, "contains", run1024};
        std::cout << myRelationship.toNode().to_json() << std::endl;

        // Create a new ID with local scope and use it to create a Record and Relationship
        axom::sina::ID myLocalID{"my_local_run", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRun{new axom::sina::Run{myLocalID, "My Sim Code", "1.2.3", "jdoe"}};
        axom::sina::Relationship myRelationship{task22, "containts", myLocalID};
    }

In the above code, the "my_local_run" ID would be replaced by a global ID on ingestion.
If this new global ID was, for example, "5Aed-BCds-23G1", then "my_local_run" would
automatically be replaced by "5Aed-BCds-23G1" in both the ``Record`` and ``Relationship``
entries.
