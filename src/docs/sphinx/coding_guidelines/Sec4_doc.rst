********************
4 Code Documentation
********************

This section contains content and formatting guidelines for the various code
documentation items mentioned in earlier sections. The aims of these 
guidelines are to:

   * Document files, data types, functions, etc. consistently.
   * Promote good documentation practices so that essential information is 
     presented clearly and lucidly, and which do not over-burden developers.
   * Generate source code documentation using the Doxygen system.


========================================
4.1 General documentation considerations
========================================

4.1.1 New source code **must** be documented following the guidelines in this section. Documentation of existing code **should** be modified to conform to these guidelines when appropriate. 

      Documentation of existing code **should** be changed when significant code
      modifications are made (i.e., beyond bug fixes and small changes) and 
      existing documentation is insufficient.

4.1.2 All header and source files **should** have comments necessary to make the code easy to understand. However, extraneous comments (e.g., documenting "the obvious") **should** be avoided.

      Code that has clear, descriptive names (functions, variables, etc.) and 
      clear logical structure is preferable to code that relies on a lot of 
      comments for understanding. To be useful, comments must be understood by 
      others and kept current with the executable code. Generally, maintenance 
      and understanding are better served by rewriting tricky, unclear code 
      than by adding comments to it.

4.1.3 End-of-line comments **should** not be used to document code logic, since they tend to be less visible than other comment forms and may be difficult to format cleanly. 

      Short end-of-line comments **may** be useful for labeling closing braces 
      associated with nested loops, conditionals, for scope in general, and 
      for documenting local variable declarations.

4.1.4 All comments, except end of line comments, **should** be indented to match the indentation of the code they are describing. Multiple line comment blocks **should** be aligned vertically on the left.

4.1.5 To make comment text clear and reduce ambiguity, code comments **should** be written in grammatically-correct complete sentences.

4.1.6 Comments **should** be clearly delimited from executable code with blank lines and "blocking characters" (see examples below) to make them stand out and, thus, improve the chances they will be read.

4.1.7 Blank lines, indentation, and vertical alignment **should** be used in comment blocks to enhance readability, emphasize important information, etc.


===================================================================
4.2 General Doxygen usage guidelines and summary of common commands
===================================================================

The Doxygen code documentation system uses C or C++ style comment sections 
with special markings and Doxygen-specific commands to extract documentation 
from source and header files. Although Doxygen provides many sophisticated 
documentation capabilities and can generate a source code manual in a variety 
of formats such as LaTeX, PDF, and HTML, these guidelines address only a small 
subset of Doxygen syntax. The goal of adhering to a simple documentation 
is that developers will be encouraged to build useful documentation when they
are writing code.

4.2.1 Doxygen comment blocks for C-only files **must** use either JavaDoc or Qt style comment block forms.

      JavaDoc style comments consist of a C-style comment block starting with 
      two \*'s, like this::

         /**
          * ...comment text...
          */

      Qt style comments add an exclamation mark (!) after the opening of a 
      C-style comment block,like this::

         /*!
          * ...comment text...
          */

      In either case, the intermediate asterisk characters ("\*" are optional, 
      but strongly encouraged.

4.2.2 Doxygen comment blocks for C++ files **must** use either JavaDoc or Qt style comment block forms (see previous item) or one of the C++ comment forms described below.

      Use a block of at least two C++ comment lines, where each line starts 
      with an additional slash::

         ///
         /// ...comment text...
         ///

      or an exclamation mark::

         //!
         //! ...comment text...
         //!

4.2.3 To be processed properly, Doxygen commands **must** be preceded with a special Doxygen character, either "\\" or "\@".

      For example, either of the following forms is acceptable Doxygen syntax 
      for providing a "brief" descriptive comment::

         \brief  ...comment text...

      or::

         @brief  ...comment text...

4.2.4 Whichever Doxygen comment block style or character used to signify a Doxygen command is used, it **must** be the same within a file.

4.2.5 Most Doxygen comments **should** appear immediately before the items they describe. 

      **Exceptions:** Inline Doxygen comments used after items such as 
      class/struct data members, enum values, function arguments, etc. **must** 
      appear after the item be **on the same line** and **must** use the 
      following syntax::

          /*!< ...comment text... */

      Note that the "<" character must appear immediately after the opening of 
      the Doxygen comment (with no space before). This tells Doxygen that the 
      comment applies to the item immediately preceding the comment. See 
      examples below.

4.2.6 A "brief" description **should** be provided in the Doxygen comment section for each of the following items: 

      * A type definition (i.e., class, struct, typedef, enum, etc.) 
      * A macro definition
      * A struct field or C++ class data member
      * A C++ class member function declaration (in the header file class 
        definition) 
      * An unbound function signature (in a header file)
      * A function implementation (when there is no description in the 
        associated header file)

      A brief comment **should** be a concise statement of purpose for an item 
      (usually no more than one line) and must start with the Doxygen command 
      "\\brief" (or "@brief").

      The Doxygen system interprets each comment as either "brief" or 
      "detailed". Brief comments appear in summary sections of the generate 
      documentation. They are typically seen before detailed comments when 
      scanning the documentation; thus good brief comments make it easier to 
      navigate a source code manual.

4.2.7 Important information of a more lengthy nature (e.g., spanning multiple lines) **should** be provided for files, major data types and definitions, functions, etc. when needed. A detailed comment **must** be separated from a brief comment with a blank line.

4.2.8 Summary of commonly used Doxygen commands

This Section provides an overview of Doxygen commands used commonly in the CS 
Toolkit source code documentation. Please see the Doxygen documentation cited
in the references at the end of these guidelines for more details and 
information about other commands that you may find useful.

Note that to be processed properly, Doxygen commands **must** be preceded with 
either "\\" or "\@" character. For brevity, we use "\\" for all commands 
described here.

   * **\\author** The "author" command (followed by a name) identifies the 
     author of a documented item. Multiple authors may be provided with each 
     on listed on its own line following the "author" keyword. 
   * **\\brief** The "brief" command is used to begin a brief description of 
     a documented item. The brief description ends at the next blank line.
   * **\\file** The "file" command is used to document a file. Doxygen requires
     that to document any global item (function, typedef, enum, etc.), the file
     in which it is defined must be documents. 
   * **\\if** and **\\endif** The "if" command, followed by a label, defines 
     the start of a conditional documentation section. The section ends with a
     matching "endif" command. Conditionals are typically used to 
     enable/disable documentation sections. For example, this may be useful if
     a project wants to provide documentation of all private class members 
     for developer documentation, but wnats to hide private members in 
     documentation for users. Conditional sections are disabled by default 
     and must be explicitly enabled in the doxygen configuration file. 
     Conditional blocks can be nested; nested sections are only enabled if 
     all enclosing sections are. The "\\elseif" command is also available to 
     provide more sophisticated control of conditional documentation.
   * **\\name** The "name" command, followed by a name containing no blank 
     spaces, is used to define a name that can be referred to elsewhere 
     in the documentation (via a link).
   * **\\param** The "param" command documents a function parameter/argument.
     It is followed by the parameter name and description. The "\\param" 
     command can be given an optional attribute to indicate usage of the 
     function argument; possible values are "[in]", "[out]", and "[in,out]".
   * **\\return** The "return" command is used to describe the return value 
     of a function.
   * **\\sa** The "sa" command (i.e., "see also") is used to refer (and 
     provide a link to) another documented item. It is followed by the target 
     of the reference (e.g., class/struct name, function name, documentation 
     page, etc.).
   * **\@{** and **\@}**  These two-character sequences begin and end a 
     grouping of documented items. Optionally, the group can be given a name 
     using the "name" command. Groups are useful for providing additional 
     organization in the documentation, and also when several items can be 
     documented with a single description (e.g., a set of simple, related 
     functions). 

   * **\\verbatim, \\endverbatim** The "verbatim/endverbatim" commands are 
     used to start/stop a block of text that is to appear exactly as it is 
     typed, without additional formatting, in the generated documentation.

   * **-** and **-#** The "-" and "-#" symbols begin an item in a bulleted 
     list or numbered list, respectively. In either case, the item ends at 
     the next blank line or next item.

   * **\\b** and **\\e** These symbols are used to make the next word bold or 
     emphasized/italicized, respectively, in the generated documentation.
   

============================
4.3 LLNL copyright statement
============================

4.3.1 Each header and source file **must** begin with a comment section containing the LLNL copyright statement (using whichever comment characters are appropriate for the programming language). For example:

.. code-block:: cpp

   /*
    * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
    * Produced at the Lawrence Livermore National Laboratory.
    *
    * All rights reserved.
    *
    * This source code cannot be distributed without permission and 
    * further review from Lawrence Livermore National Laboratory.
    */


============================
4.4 File documentation
============================

4.3.1 Each header files that declares unbound functions, defines enums, typedefs, etc. **must** have a Doxygen file prologue similar to the following:

.. code-block:: cpp

   /*!
    ***************************************************************************
    *
    * \file ...optional name of file...
    *
    * \brief A brief statement describing the file contents/purpose. (optional)
    *
    * Optional detailed explanatory notes about the file.
    *
    * \author Name of file author (optional)
    *
    ****************************************************************************
    */

4.3.2 The Doxygen command "\\file" **must** appear first in the file prologue.

      The "\\file" command identifies the comment section as documentation 
      for the file. Doxygen requires that the file itself must be documented 
      for documentation to be generated for any global item (global function, 
      typedef, enum, etc.) defined in the file.

      The file name may include (part of) the path if the file name is not 
      unique. If the file name is omitted on the line after the "\\file" 
      command, then any documentation in the comment block will belong to 
      the file in which it is located instead of the summary documentation 
      in the listing of documented files.

4.3.3 A brief statement of purpose for the file **should** appear as the first comment after the file. If included, the brief statement, **must** be preceded by the "\\brief" command.

      Brief documentation statements are often helpful to those scanning the 
      documentation.

4.3.4 Any detailed notes about the file **may** be included after the brief comment. If this is done, the detailed comments **must** be separated from the brief statement by a blank line.

4.3.4 The name of the original author of the file **may** be entered after the file notes. If the author's name is included, it **must** be preceded by the "\\author" command.


========================
4.5 Type documentation
========================

4.5.1 Each type definition (i.e., class, struct, enum, typedef, etc.) and macro definition appearing in a header file **must** have a Doxygen type definition comment prologue immediately before it. For example

.. code-block:: cpp

   /*!
    ****************************************************************************
    *
    * \brief A brief statement of purpose of the type or macro.
    *
    * Optional detailed information that is helpful in understanding the
    * purpose, usage, etc. of the type/macro ...
    *
    * \sa optional cross-reference to other types, functions, etc...
    * \sa etc...
    *
    * \warning This class is only partially functional.
    *
    ****************************************************************************
    */

Note that Doxygen requires that a compound entity, such as a class, struct, 
etc. must be documented in order to document any of its members.

4.5.2 A brief statement describing the type **must** appear as the first text comment using the Doxygen command "\\brief".

4.5.3 Important details about the item **should** be included after the brief comment and, if included, **must** be separated from the brief comment by a blank line.

4.5.4 Cross-references to other items, such as relevant major types, important functions, etc., **should** be included at the end of the prologue to enhance the navigability of the Doxygen documentation. 

      The Doxygen command "\\sa" (for "see also") **should** appear before each
      such cross-reference so that links are generated in the documentation.

4.5.6 Caveats or limitations about the documented type **should** be noted using the "\\warning" Doxygen command as shown above.


===============================
4.6 Data member documentation
===============================

4.6.1 Each struct field, C++ class data member, etc. **should** have a descriptive comment indicating its purpose. 

     This comment may as appear as a prologue before the item, such as::

        /*!
         *
         * \brief Brief statement of purpose of data member m_mode.
         *
         * Optional detailed information about m_mode...
         */
        int m_mode;

     or, it may appear after the item as an inline comment such as::

        int m_mode; /*!< \brief Brief statement of purpose of m_mode... */

4.6.2 Regardless of which documentation form is used, a brief description of purpose of the definition **must** be included using the Doxygen command "\\brief".

4.6.3 When documenting a data item inline (as in the second example above), the comment must follow the item on the same line.

     The form of an inline Doxygen comment is::

         /*!< \brief ...comment text... */

     Note that the "<" character must be included immediately after the start 
     of the Doxygen comment form (with no space between). This tells Doxygen 
     that the comment corresponds to the item immediately preceding it.

4.6.4 Any detailed notes about an item, if included, **must** appear after the brief comment and be separated from the brief comment with a blank line. 

4.6.5 When a detailed comment is provided, or the brief statement requires more than one line, the prologue comment form **should** be used instead of the inline form to make the documentation easier to read.

4.6.6 If the names of data members are sufficiently clear that their meaning and purpose are obvious to other developers (which should be determined in a code review), then the members **may** be grouped together and documented with a single descriptive comment.

      An example of Doxygen syntax for such a grouping is::

         //@{
         //!  @name Data member description...

         int m_member1;
         int m_member2;
         ...
         //@}


==========================
4.7 Function documentation
==========================

4.7.1 Each unbound functions **should** be be documented with a function prologue in the header file where its prototype appears or in a source file immediately preceding its implementation.

4.7.2 Since C++ class member member functions define the class interface, they **should** be documented with a function prologue immediately preceding their declaration in the class definition.

The following examples show two function prologue variations that may 
be used to document a method in a class definition. The first shows how
to document the function arguments in the function prologue.

.. code-block:: cpp

      /*!
       *************************************************************************
       *
       * \brief Initialize a Foo object to given operation mode.
       *
       * The "read" mode means one thing, while "write" mode means another.
       *
       * \return bool indicating success or failure of initialization.
       *              Success returns true, failure returns false.
       *
       * \param[in] mode OpMode enum value specifying initialization mode.
       *                 ReadMode and WriteMode are valid options.
       *                 Any other value generates a warning message and the
       *                 failure value ("false") is returned.
       *
       *************************************************************************
       */
       bool initMode(OpMode mode);

The second example shows how to document the function arguments inline.

.. code-block:: cpp

      /*!
       ************************************************************************
       *
       * @brief Initialize a Foo object to given operation mode.
       *
       * The "read" mode means one thing, while "write" mode means another.
       *
       * @return bool value indicating success or failure of initialization.
       *             Success returns true, failure returns false.
       *
       *************************************************************************
       */
       bool initMode(OpMode mode /*!< [in] ReadMode, WriteMode are valid options */ );

Note that the first example uses the "\\" character to identify Doxygen 
commands; the second uses "@". Also, the "<" character must appear immediately 
after the start of the Doxygen comment form (with no space between). This 
tells Doxygen that the comment corresponds to the item immediately preceding it.

4.7.3 A brief statement of purpose for a function must appear as the first text comment after the Doxygen command "\\brief" (or "@brief"). 

4.7.4 Any detailed notes about a function, when included, **must** appear after the brief comment and **must** be separated from the brief comment by a blank line.

4.7.4 If the function has a non-void return type, the return value **should** be documented in the prologue using the Doxygen command "\return" (or "@return") preceding a description of the return value. 

      Functions with "void" return type and C++ class constructors and 
      destructors **should not** have such documentation.

4.7.5 Function arguments **should** be documented in the function prologue or inline (as shown above) when the intent or usage of the arguments is not obvious. 

      The inline form of the comment may be preferable when the argument 
      documentation is short. When a longer description is provided (such as 
      when noting the range of valid values, error conditions, etc.) the 
      description **should** be placed within the function prologue for 
      readability. However, the two alternatives for documenting function 
      arguments **must not** be mixed within the documentation of a single 
      function to reduce confusion. 

      In any case, superfluous documentation should be avoided. For example, 
      when there are one or two arguments and their meaning is obvious from 
      their names or the description of the function, providing no comments is 
      better than cluttering the code by documenting the obvious. Comments 
      that impart no useful information are distracting and less useful than 
      no comment at all.

4.7.6 When a function argument is documented in the prologue comment section, the Doxygen command "\param" **should** appear before the comment as in the first example above.

4.7.7. The "in/out" status of each function argument **should** be documented.

       The Doxygen "\param" command supports this directly by allowing such an
       attribute to be specified as "\param[in]", "\param[out]", or 
       "\param[in,out]". Although the inline comment form does not support 
       this, such a description **should** be included; e.g., by using "[in]", 
       "[out]", or "[in,out]" in the comment.

4.7.8 Short, simple functions (e.g., inline methods) **may** be grouped together and documented with a single descriptive comment when this is sufficient.

      An example of Doxygen syntax for such a grouping is::

         //@{
         //! @name Setters for data members

         void setMember1(int arg1) { m_member1 = arg1; }
         void setMember2(int arg2) { m_member2 = arg2; }

         //@}

4.7.9 Typically, important implementation details about a function **should** be documented in the source file where the function is implemented. 

      Header file documentation **should** include only purpose and usage 
      information germane to an interface. When a function has separate 
      implementation documentation, the comments **must** not contain Doxygen 
      syntax. Using Doygen syntax to document an item in more than one location 
      (e.g., header file and source file) can cause undesired Doxygen 
      formatting issues and potentially confusing documentation.
      

      A member of a class may be documented as follows in the source file 
      for the class as follows::

        /*
         ***********************************************************************
         *
         * Set operation mode for a Foo object.
         *
         * Important detailed information about what the function does...
         *
         ***********************************************************************
         */
         bool Foo::initMode(OpMode mode)
         {
            ...function body...
         }
