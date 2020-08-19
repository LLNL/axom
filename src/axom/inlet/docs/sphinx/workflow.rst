.. _inlet_workflow_label:

Workflow
========

Inlet's workflow is broken into the three following steps.

.. _inlet_defining_schema_label:

Defining Schema
---------------

The first step in using Inlet is to defining the schema of your input deck.
Inlet defines an input deck into two basic classes: Tables and Fields. Basically
Fields are individual values and Tables hold groups of Fields and Tables.

Define the schema by using the following functions, on either the main Inlet class, for
global Tables and Fields, or on individual Table classes, for Tables and Fields under that Table:

========================= ===================
Name                      Description
========================= ===================
addTable                  Adds a Table to the input deck schema with the given name.
addBool                   Adds a boolean Field to the global or parent Table with the given name.
addDouble                 Adds a double Field to the global or parent Table with the given name.
addInt                    Adds a integer Field to the global or parent Table with the given name.
addString                 Adds a string Field to the global or parent Table with the given name.
========================= ===================

All possible Tables and Fields that are can be found in the input deck should be defined
at this step.  Use the ``required`` class member function on the Table and Field class to indicate that
they have to present in the given input deck. For example, the following indicates that
the given input deck must defined a global Field named ``dimensions``:

.. code-block:: c

    myInlet->addInt("dimensions")->required(true);


The value of the Field is read and stored into the Sidre datastore when you call the appropriate
add function. You can also set a default value to each field via the type-safe ``Field::defaultValue()``
member functions. Doing so will populate the corresponding Fields value if the specific Field is not
present in the input deck. For example, the follwing will set the ``foo`` field value to true if it is 
not present in the given input deck:

.. code-block:: c

    myInlet->addBool("foo")->defaultValue(true);


.. _inlet_verification_label:

Verification
------------

This step helps ensure that the given input deck follows the rules expected by the code.  These
rules are not verified until you call ``Inlet::verify()``.  Doing so will return true/false and
output SLIC warnings to indicate which Field or Table violated which rule.

As shown above, both Tables and Fields can be marked as ``required``. Fields have two additional
basic rules that can be enforced with the following ``Field`` class member functions:

========================= ===================
Name                      Description
========================= ===================
validValues               Indicates the Field can only be set to a valid values.
range                     Indicates the Field can only be set to inclusively between two values.
========================= ===================

Inlet also provides functionality to write your own custom rules via callable lambda verifiers.
Fields and Tables can both register one lambda each via their ``registerVerifier()`` member functions.
The following example adds a custom verifier that simply verifies that the dimension of the simulation
match up with the dimensions of a given vector:

.. code-block:: c

    void example() {
    auto lr = std::make_shared<axom::inlet::LuaReader>();
    lr->parseString("dimensions = 2; vector = { x = 1; y = 2; z = 3; }");
    axom::sidre::DataStore ds;
    auto myInlet = std::make_shared<axom::inlet::Inlet>(lr, ds.getRoot());
    myInlet->addInt("dimensions")->required(true);
    auto v = myInlet->addTable("vector")->required(true);

    v->addInt("x");

    v->registerVerifier([&]() -> bool {
        int dim;
        myInlet->get("dimensions", dim);
        int value;  // field value doesnt matter just that it is present in input deck
        bool x_present = v->hasChildField("x") && myInlet->get("vector/x", value);
        bool y_present = v->hasChildField("y") && myInlet->get("vector/y", value);
        bool z_present = v->hasChildField("z") && myInlet->get("vector/z", value);
        if(dim == 1 && x_present) {
        return true;
        }
        else if(dim == 2 && x_present && y_present) {
        return true;
        }
        else if(dim == 3 && x_present && y_present && z_present) {
        return true;
        }
        return false;
    });

    // We expect verification to be unsuccessful since the only Field
    // in vector is x but 2 dimensions are expected
    myInlet->verify() ? std::cout << "Verification was successful\n" 
                        : std::cout << "Verification was unsuccessful\n";

    v->addInt("y");
    std::cout << "After adding a new dimension:\n";

    // We expect the verification to succeed because vector now contains
    // both x and y to match the 2 dimensions
    myInlet->verify() ? std::cout << "Verification was successful\n" 
                        : std::cout << "Verification was unsuccessful\n";

    }

.. note::  ``Inlet::getGlobalTable()->registerVerifier()`` can be used to add a verifier to apply rules
  to the Fields at the global level.


.. _inlet_accessing_data_label:

Accessing Data
--------------

Accessing data happens with some functions.. docs incoming.
