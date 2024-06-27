.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _documents-label:

=========
Documents
=========

Sina ``Document`` objects are a way to represent the top-level object of a
JSON file that conforms to the Sina schema. When serialized, these documents
can be ingested into a Sina database and used with the Sina tool.

``Document`` objects follow a very basic JSON layout consisting of two entries:
"records" and "relationships". Each of these entries will store a list of their
respective objects. An example of an empty document is shown below:

.. code:: json

    {
        "records": [],
        "relationships": []
    }

The "records" list can contain ``Record`` objects and their inheritying types,
such as ``Run`` objects. The "relationships" list can contain ``Relationship``
objects. For more information on these objects, see `Records <./records>`_
and `Relationships <./relationships>`_.

--------------------
Assembling Documents
--------------------

``Document`` objects can be assembled programatically. To accomplish this:

1. Create a new instance of the ``Document`` class
2. Create a ``Record``
3. Add the instance of the ``Record`` with the ``add`` method

On the `Sina C++ User Guide <./index>`_ page, you can see an example of this
process. Below we will expand on this example to add a ``Relationship``:

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        // Create a new document
        axom::sina::Document document;

        // Create a record of this specific study
        // This study has an ID of "study1", which has to be unique to this file
        axom::sina::ID studyID{"study1", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> study{
            new axom::sina::Record{studyID, "UQ study"}};

        // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
        // The run has an ID of "run1", which has to be unique to this file.
        axom::sina::ID runID{"run1", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> run{
            new axom::sina::Run{runID, "My Sim Code", "1.2.3", "jdoe"}};

        // Create a relationship between the study and the run
        // Here we're saying that the study contains the run
        axom::sina::Relationship relationship{studyID, "contains", runID};

        // Add the run, study record, and relationship to the document
        document.add(std::move(run));
        document.add(std::move(study));
        document.add(relationship);

        // Save the document directly to a file.
        saveDocument(document, "MySinaData.json");
    }

After executing the above code, the resulting ``MySinaData.json`` file will
look like so:

.. code:: json

    {
        "records": [
            {
                "type": "run",
                "local_id": "run1",
                "application": "My Sim Code",
                "version": "1.2.3",
                "user": "jdoe"
            },
            {
                "type": "UQ study",
                "local_id": "study1"
            }
        ],
        "relationships": [
            {
                "predicate": "contains",
                "local_subject": "study1",
                "local_object": "run1"
            }
        ]
    }

------------------------------
Generating Documents From JSON
------------------------------

Alternatively to assembling ``Document`` instances programatically, it is
also possible to generate ``Document`` objects from existing JSON files
or JSON strings.

Using our same example from the previous section, if we instead had the
``MySinaData.json`` file prior to executing our code, we could generate
the document using Sina's ``loadDocument()`` function:

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        axom::sina::Document myDocument = axom::sina::loadDocument("MySinaData.json");
    }

Similarly, if we had JSON in string format we could also load an instance
of the ``Document`` that way:

.. code:: cpp

    #include "axom/sina.hpp"

    int main (void) {
        std::string my_json = "{\"records\":[{\"type\":\"run\",\"id\":\"test\"}],\"relationships\":[]}";
        axom::sina::Document myDocument = axom::sina::Document(my_json, axom::sina::createRecordLoaderWithAllKnownTypes());
        std::cout << myDocument.toJson() << std::endl;
    }

---------------------------------------------------------
Obtaining Records & Relationships from Existing Documents
---------------------------------------------------------

Sina provides an easy way to query for both ``Record`` and ``Relationship``
objects that are associated with a ``Document``. The ``getRecords()`` and
``getRelationships()`` methods will handle this respectively.

Below is an example showcasing their usage:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create a new document
        axom::sina::Document document;

        // Create a record of this specific study
        // This study has an ID of "study1", which has to be unique to this file
        axom::sina::ID studyID{"study1", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> study{
            new axom::sina::Record{studyID, "UQ study"}};

        // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
        // The run has an ID of "run1", which has to be unique to this file.
        axom::sina::ID runID{"run1", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> run{
            new axom::sina::Run{runID, "My Sim Code", "1.2.3", "jdoe"}};

        // Create a relationship between the study and the run
        // Here we're saying that the study contains the run
        axom::sina::Relationship relationship{studyID, "contains", runID};

        // Add the run, study record, and relationship to the document
        document.add(std::move(run));
        document.add(std::move(study));
        document.add(relationship);

        // Query for a list of records and relationships
        auto &records = document.getRecords();
        auto &relationships = document.getRelationships();

        std::cout << "Number of Records: " << records.size() << std::endl;
        std::cout << "Number of Relationships: " << relationships.size() << std::endl;
    }

Running this will show that both records and the one relationship were
properly queried:

.. code:: bash

    Number of Records: 2
    Number of Relationships: 1
