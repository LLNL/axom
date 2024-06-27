.. ## Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
.. ## other Axom Project Developers. See the top-level LICENSE file for details.
.. ##
.. ## SPDX-License-Identifier: (BSD-3-Clause)

.. _records-label:

=======
Records
=======

Sina ``Record`` objects are used to represent the data to be stored for your
study. Some examples for the natural scope of ``Record`` objects include things
like:

    - a single run of an application
    - an msub job
    - a cluster of runs that has some metadata attached to the cluster (this
      ``Record`` might have a "contains" :doc:`Relationship <relationships>` for all
      the runs within it)

Each ``Record`` must have an ID and a type that you define. Additionally, Each
``Record`` can have a list of ``File`` objects and a list of ``Datum`` objects.

.. contents:: Topics Covered in this Page
   :depth: 2
   :local:

------
Datums
------

Sina ``Datum`` objects help track the value and (optionally) tags and/or units of a
value associated with a ``Record``. In the Sina schema, a ``Datum`` always belongs
to a ``Record`` or one of ``Record``'s inheriting types.

Some examples of potential ``Datum`` values would be:

    - a scalar
    - a piece of metadata
    - an input parameter

The value of a ``Datum`` may be a string, a double, an array of strings, or an array
of doubles.

Below showcases an example of creating an instance of ``Datum`` with an array of
strings and adding it to a ``Record``:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create the record
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Create the datum with an array of strings
        std::vector<std::string> myTags{"input"};
        axom::sina::Datum myDatum{12, myTags};

        // Add the datum to the record
        myRecord->add("my_scalar", std::move(myDatum));
        std::cout << myRecord->toNode().to_json() << std::endl;
    }

Once executed, this code will output:

.. code:: json

    {
        "local_id": "my_record",
        "type": "my_type",
        "data":
        {
            "my_scalar":
            {
                "tags": ["input"],
                "value":12.0
            }
        }
    }

.. _datum-type-label:

+++++++++++++++++++
Checking Datum Type
+++++++++++++++++++

It's possible to check the type of a ``Datum`` with the ``getType()`` method. Types
are tracked in an enumeration called ``ValueType``. The enumeration is as follows:

    - 0: string
    - 1: scalar
    - 2: array of strings
    - 3: array of scalars

Below is an example of this in action:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Define 3 different datums
        axom::sina::Datum myDatum{12.34};
        std::string value = "foobar";
        axom::sina::Datum myOtherDatum{value};
        std::vector<double> scalars = {1, 2, 20.0};
        axom::sina::Datum myArrayDatum{scalars};

        // Prints 0, corresponding to string
        std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(myDatum.getType()) << std::endl;

        // Prints 1, corresponding to scalar
        std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(myOtherDatum.getType()) << std::endl;

        // Prints 3, corresponding to scalar array
        std::cout << static_cast<std::underlying_type<axom::sina::ValueType>::type>(myArrayDatum.getType()) << std::endl;
    }

++++++++++++++++++++++
Setting Units and Tags
++++++++++++++++++++++

For certain ``Datum`` instances it may be helpful to assign them units and/or tags.
This can be accomplished with the ``setUnits()`` and ``setTags()`` methods respectively.

Below is an example of this functionality:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Define 2 different datums
        axom::sina::Datum myDatum{12.34};
        std::vector<double> scalars = {1, 2, 20.0};
        axom::sina::Datum myArrayDatum{scalars};

        // Set the units for one datum and the tags for the other
        myDatum.setUnits("km/s");
        std::vector<std::string> tags = {"input", "core"};
        myArrayDatum.setTags(tags);
    }

+++++++++++++++++++++++++++++++++++++
Viewing Datum From an existing Record
+++++++++++++++++++++++++++++++++++++

Sometimes it's necessary to obtain the current ``Datum`` instances from an existing
``Record``. To do this, you can utilize the ``Record`` object's ``getData`` method.
This method will return an unordered map of ``Datum`` instances.

Below is an example of this process:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Define 3 different datums
        axom::sina::Datum myDatum{12.34};
        std::string value = "foobar";
        axom::sina::Datum myOtherDatum{value};
        std::vector<double> scalars = {1, 2, 20.0};
        axom::sina::Datum myArrayDatum{scalars};

        // Create a record to store the datum
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Add the datum instances to the record
        myRecord->add("datum1", std::move(myDatum));
        myRecord->add("datum2", std::move(myOtherDatum));
        myRecord->add("datum3", std::move(myArrayDatum));

        // Query the datum
        auto &data = myRecord->getData();

        // Print the keys and type of datum
        for (const auto& pair : data) {
            std::cout << pair.first << " is type: " << static_cast<std::underlying_type<axom::sina::ValueType>::type>(pair.second.getType()) << std::endl;
        }
    }

Executing this code will print out:

.. code:: bash

    datum1 is type: 1
    datum2 is type: 0
    datum3 is type: 3

Which, we know from `Checking Datum Type <#datum-type-label>`_, signifies that
datum1 is a scalar, datum2 is a string, and datum3 is an array of scalars.

Using this knowledge we can modify our code to show us the current datum values:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Define 3 different datums
        axom::sina::Datum myDatum{12.34};
        std::string value = "foobar";
        axom::sina::Datum myOtherDatum{value};
        std::vector<double> scalars = {1, 2, 20.0};
        axom::sina::Datum myArrayDatum{scalars};

        // Create a record to store the datum
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Add the datum instances to the record
        myRecord->add("datum1", std::move(myDatum));
        myRecord->add("datum2", std::move(myOtherDatum));
        myRecord->add("datum3", std::move(myArrayDatum));

        // Query the datum
        auto &data = myRecord->getData();

        // Print the datum values
        std::cout << "datum1: " << data.at("datum1").getScalar() << std::endl;
        std::cout << "datum2: " << data.at("datum2").getValue() << std::endl;
        std::cout << "datum3: ";
        for (const auto& value : data.at("datum3").getScalarArray()) {
            std::cout << value << " ";
        }
        std::cout << std::endl;
    }

This will provide the following output:

.. code:: bash

    datum1: 12.34
    datum2: foobar
    datum3: 1 2 20 

-----
Files
-----

Sina ``File`` objects help track the location (URI) and mimetype of a file on the
file system, plus any tags. In the Sina schema, a ``File`` always belongs to a ``Record``
or one of ``Record``'s inheriting types.

Every ``File`` must have a URI, while mimetype and tags are optional.

Below is an example showcasing how to create a file and add it to a record:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create 2 different files
        axom::sina::File myFile{"/path/to/file.png"};
        myFile.setMimeType("image/png");
        axom::sina::File myOtherFile{"/path/to/other/file.txt"};
        myOtherFile.setTags({"these", "are", "tags"});

        // Create a record to store the files
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Add the files to the record
        myRecord->add(myFile);
        myRecord->add(myOtherFile);

        std::cout << myRecord->toNode().to_json() << std::endl;
    }

This code will produce the following output:

.. code:: json

    {
        "type": "my_type",
        "local_id": "my_record",
        "files": 
        {
            "/path/to/other/file.txt": 
            {
                "tags": 
                [
                    "these",
                    "are",
                    "tags"
                ]
            },
            "/path/to/file.png":
            {
                "mimetype": "image/png"
            }
        }
    }

Similarly, files can be removed from a ``Record`` with the ``remove()`` method:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create 2 different files
        axom::sina::File myFile{"/path/to/file.png"};
        myFile.setMimeType("image/png");
        axom::sina::File myOtherFile{"/path/to/other/file.txt"};
        myOtherFile.setTags({"these", "are", "tags"});

        // Create a record to store the files
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Add the files to the record
        myRecord->add(myFile);
        myRecord->add(myOtherFile);

        // Remove a file from the record
        myRecord->remove(myFile);

        std::cout << std::cout << myRecord->toNode().to_json() << std::endl;
    }

As we see from the output, the contents of ``myFile`` are no longer in the
``Record`` instance:

.. code:: json

    {
        "type": "my_type",
        "local_id": "my_record",
        "files": 
        {
            "/path/to/other/file.txt": 
            {
                "tags": 
                [
                    "these",
                    "are",
                    "tags"
                ]
            }
        }
    }

+++++++++++++++++++++++++++++++++++++
Viewing Files From an Existing Record
+++++++++++++++++++++++++++++++++++++

Sometimes it's necessary to view the current ``File`` instances that are stored in
a ``Record``. This can be accomplished by using the ``Record`` object's ``getFiles()``
method which returns an unordered map of ``File`` instances.

Below is an expansion of the previous example where we query the ``Record`` instance
for files:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create 2 different files
        axom::sina::File myFile{"/path/to/file.png"};
        myFile.setMimeType("image/png");
        axom::sina::File myOtherFile{"/path/to/other/file.txt"};
        myOtherFile.setTags({"these", "are", "tags"});

        // Create a record to store the files
        axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
        std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

        // Add the files to the record
        myRecord->add(myFile);
        myRecord->add(myOtherFile);

        // Query the record for files
        auto &files = myRecord->getFiles();
        for (const auto& file : files) {
            std::cout << "File with URI '" << file.getUri() << "' has mimetype '" << file.getMimeType() << "' and tags '";
            for (const auto& tag : file.getTags()) {
                std::cout << tag << " ";
            }
            std::cout << "'" << std::endl;
        }
    }

The above code will output:

.. code:: bash

    File with URI '/path/to/file.png' has mimetype 'image/png' and tags ''
    File with URI '/path/to/other/file.txt' has mimetype '' and tags 'these are tags '

----
Runs
----

Sina ``Run`` objects are a subtype of Sina ``Record`` objects corresponding to
a single run of an application as specified in the Sina schema. Similar to ``Record``
objects, ``Run`` objects also require an ID; however, they do *not* require a type as
it will automatically be set to "run". In addition to IDs, there are a few other
required fields:

    - application: the application/code used to create the ``Run``
    - version: the version of the application used to create the ``Run``
    - user: the username of the person who ran the application that generated this ``Run``

The below code shows an example of how a ``Run`` can be created:

.. code:: cpp

    #include "axom/sina.hpp"

    int main(void) {
        // Create the ID first
        axom::sina::ID run1ID{"run1", axom::sina::IDType::Local};

        // Create a run with:
        // ID: run1ID
        // application: "My Sim Code"
        // version: "1.2.3"
        // user: "jdoe"
        std::unique_ptr<sina::Record> run1{new sina::Run{run1ID, "My Sim Code", "1.2.3", "jdoe"}};
    }
