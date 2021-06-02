.. _inlet_verification_page_label:

############
Verification
############

Before input file data can be accessed, it must first be verified by calling the ``verify()``
method of the top-level ``Inlet`` object. This will return a ``bool`` indicating whether the
provided input conformed to the schema and specific violations of the schema are logged via
SLIC by default.  If you would like to suppress the SLIC warnings and process the list of
verification errors instead, you can pass a ``std::vector<inlet::VerificationError>`` to the
``verify()`` method as follows:

.. code-block:: C++

  std::vector<inlet::VerificationError> errors;
  inlet.verify(&errors);

You can then iterate over the list of errors, each of which contains the path within the input file
of the offending ``Container``, ``Field``, or ``Function`` and the corresponding message.

This section describes the verification rules that apply to each possible element of the Inlet
hierarchy, namely, ``Container``, ``Field``, and ``Function``.

Container Verification
----------------------

If a ``Container`` is marked as required (via the ``required()``) method, then the ``Container`` must
have a ``Field`` or ``Function`` that was present in the input or contain a sub-``Container`` that does.
This does not apply to a ``Container`` that corresponds to an array or dictionary, as empty collections
provided by the user are considered valid.  Consider the following definition and input file:

.. code-block:: C++

  addIntArray("foo").required();

.. code-block:: Lua

  foo = { }

Inlet verification will succeed for the above snippets.

If a ``Container`` corresponds to an array or dictionary, the elements of the array must all be of the requested
type, if any were provided.  This restriction applies even if the array/dictionary was not marked as ``required``.

If a verification function was provided via the ``registerVerifier()`` method, this function must 
return ``true`` when passed the corresponding ``Container`` object.

.. note::
  Since a ``Container`` represents internal nodes in the Inlet hierarchy, its verification status is
  dependent on that of its child objects.  A ``Container`` is only considered valid if all of its child
  ``Container``, ``Field``, and ``Function`` objects are also valid.

Field Verification
------------------

If a ``Field`` is marked as required (via the ``required()``) method, then a value must be provided in the input.

Providing a value of the wrong type will result in verification failure, even if the field was not marked as ``required``.

If a range (inclusive upper/lower bounds) of valid values is specified with the ``range()`` method, both the provided value
(if applicable) and default value (if specified with ``defaultValue()``) must fall within the range.

If a set of valid values is specified with the ``validValues()`` method, both the provided value
(if applicable) and default value (if specified with ``defaultValue()``) must be included in the set.

If a verification function was provided via the ``registerVerifier()`` method, this function must 
return ``true`` when passed the corresponding ``Field`` object.

Function Verification
---------------------

If a ``Function`` is marked as required (via the ``required()``) method, then a function must be provided in the input.

If a verification function was provided via the ``registerVerifier()`` method, this function must 
return ``true`` when passed the corresponding ``Function`` object.

Unexpected Entries in Input Files
---------------------------------

In order to better detect user error, e.g., misspelled names, Inlet provides a method to retrieve the names of entries
in the input file that were not requested in the schema definition phase.  Given a top-level ``Inlet`` object named ``inlet``,
the names can be retrieved as follows:

.. code-block:: C++

  std::unordered_set<std::string> unexpected_names = inlet.unexpectedNames();
