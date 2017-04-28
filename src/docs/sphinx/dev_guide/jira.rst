.. ##
.. ## Copyright (c) 2016, Lawrence Livermore National Security, LLC.
.. ##
.. ## Produced at the Lawrence Livermore National Laboratory.
.. ##
.. ## All rights reserved.
.. ##
.. ## This file cannot be distributed without permission and
.. ## further review from Lawrence Livermore National Laboratory.
.. ##

.. _issuetracking-label:

======================================================
Issue Tracking and Work Planning (JIRA)
======================================================


.. _releasecycle-label:

----------------------------------------------
Release Cycles and Development Planning
----------------------------------------------

We use our `JIRA project space <https://lc.llnl.gov/jira/browse/ATK>`_ for
issue tracking and project planning. In JIRA, you can create issues, edit
them, comment on them, check issue status, group them together for sprint
development, and search for issues in various ways, including setting up
filters to customize your searches.


.. note:: Fill this in....


.. _issueworkflow-label:

------------------------
Issue Workflow (JIRA)
------------------------

Issue states
^^^^^^^^^^^^^

We have customized our issue workflow to make it simple and easy to understand.
Specifically, each issue has three possible states:

  Open.
    Every issues starts out as an open issue. An open issue can
    be assigned to someone or unassigned. When an issue is assigned, this
    means that the assignee owns the issue and is responsible for working
    on it. An open issue that is unassigned has not been been discussed or
    reviewed, or no decision to act on it has been made. Typically, an open
    issue means that it is not being worked on.
  In Progress.
    An issue in progress is one that is actively being worked on.
  Closed.
    When an issue is closed, work on it has been completed, or
    a decision has been made that it will not be addressed.

An open issue can transition to either in progress (work has started on it)
or closed. An in progress issue can transition to either open (work on it
has stopped, but it is not finished) or closed. Finally, a closed issue
can be re-opened, which changes its state to open.

The figure below shows issue state transitions in our JIRA workflow.

.. figure:: jira-issue.png

   This figure shows allowed state transitions in our JIRA issue workflow.


Creating a new issue
^^^^^^^^^^^^^^^^^^^^^

To create a new issue, click the 'Create' button at the top of the Axom
JIRA project page and enter information in the issue fields. Filling in the
fields properly greatly helps other team members search through project issues
to find what they are looking for. Note that issue fields marked with a red
asterisk are required. The others are not required, but may be used to include
helpful information. The main issues we use regularly are:

  Project
    Axom will show up as the default. You shouldn't need
    to change this.
  Issue Type
    We use only three issue types: *Bug*, *New Feature*, and
    *Task*. A bug is something broken that needs to be fixed. A new feature
    is something to add that increases functionality, enhances an interface,
    etc. Task is a "catch-all" issue type for any other issue.
  Summary
    Provide a short descriptive summary. A good (and brief)
    summary makes it easy to scan a list of issues to find one you are
    looking for.
  Priority
    Select an appropriate issue priority to impart its level
    of importance or urgency. Clicking on the question mark to the right of
    the priority field provides a description of each option.
  Components
    Each issue is labeled with the Axom component it
    applies to. Other "component" labels indicate build system issues,
    documentation issues, etc.
  Assignee
    Unless you are certain which team member should be assigned
    the issue, choose 'Unassigned'. This will indicate that the issue requires
    discussion and review before it is assigned. The default assignee is the
    owner of the component you chose earlier if you make no choice.
  Reporter
    Unless you explicitly enter someone in this field, you, as
    the issue creator, will be the reporter. This is the correct choice in
    almost all cases.
  Description
    The description field should be used to include important
    details about the issue that will help the developer who will work on it.
  Environment
    The environment field can be useful when an issue affects a particular
    compiler or platform.

You may also use the other fields that appear if you think they will help
describe the issue. However, the team seldom uses fields apart from the list
above.

Starting and stopping work on an issue
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When you begin work on an issue, you should open it, make sure it is
assigned to you and click the 'Start Progress' button at the top of the issue.
This changes its status to *In progress*.

If there is still work to do on the issue, but you will stop working on it
for a while, you can click the 'Stop Progress' button at the top of the
issue. This changes its status back to *Open*.

Closing an issue
^^^^^^^^^^^^^^^^^

When work is completed on an issue (which includes testing, adding
new documentation if needed, etc.), or the issue will not be addressed,
it should be closed. To close an issue, click the 'Close' button and select
the appropriate issue resolution. There are two options: *Done* and *Won't Fix*.
'Done' means that the issue is resolved. 'Won't Fix' means that the issue will
not be addressed for some reason.

When closing an issue, adding information to the 'Comment' field may be
helpful. For example, when an issue is closed as 'Won't Fix', it is helpful to
enter a brief explanation as to why this is so.

Issue assignee
^^^^^^^^^^^^^^^

Note that an assigned issue can be assigned to someone else to work on it.
An assigned issue can also be set back to 'Unassigned' if it needs further
discussion by the team.

JIRA tips
^^^^^^^^^^

Here are some links to short videos (a couple of minutes each) that
demonstrate how to use JIRA features:

   * `JIRA Instant Search Bar Demo <https://www.youtube.com/watch?v=ZmACxhzXLco&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=3>`_
   * `JIRA System Files Demo <https://www.youtube.com/watch?v=O08oySq043w&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=4>`_
   * `Creating and Editing JIRA Issues <https://www.youtube.com/watch?v=EsQ__dR6Nrw&list=PLlALqRAjvdnGB_T0GAB1Fk2rVZgnJJAOa&index=5>`_
