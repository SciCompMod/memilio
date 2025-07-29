Developer workflow
========================

We are always happy about contributions to the project! Here you can find more information on our coding guidelines, our git workflow, benchmarking our models and writing documentation. 

Coding guidelines
---------------------

All software is built in modules, unit tests have to be added for each module/functionality.

The CI pipeline also automates some code style enforcement via a ``pre-commit``.
We recommend to configure it locally such that it runs automatically on every commit:

.. code:: bash

    pip install pre-commit
    pre-commit install


For more information about ``pre-commit`` check `here <https://docs.pymc.io/en/latest/contributing/python_style.html>`_ and this short video series: https://calmcode.io/pre-commit/the-problem.html

Please be aware that the ``isort`` pre-commit hook accidentally sorts our own code with third party libraries, also see: https://github.com/PyCQA/isort/issues/2068 . Be therefore sure to not commit python code from a worktree.

C++ Coding guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



C++ Standard:

 - C++ 20, other standards are currently not supported.

Namespaces:

  - Use the ``mio`` namespace. 

Naming rules:

  - Classes begin with large Letters , e.g. ``class MyClass``
  - functions, methods, variables use small letters + underscore, e.g. ``my_awesome_function`` 
  - member variables should be generally private (we allow exceptions from this rule) and should be named with a leading ``m_``, e.g. ``m_my_member``.

Return Values:

  - If only one object is output, use return, for multiple objects, pass by reference (we still have to check ``std::expected``)
  - The semantics of return value arguments have to make clear, how the ownership is handled

    - If the function creates an object (allocates), pass it as ``std::unique_ptr<T>&``
    - If the function simply changes an object, pass is as ``T&``

  - Avoid producing unnecessarily long outputs. Ensure that all output is concise and limited to relevant information.

Exceptions:

  - In order to avoid MPI deadlocks, do not use exceptions. Use logging or return codes.

Logging:

  - Do not use printfs
  - Use the logging functions from ``logging.h``
  - For debug logs, use ``mio::log_debug(msg)``

Includes:

  - Please use include guards with capitalized name of the header file (``test.h -> #ifndefine TEST_H``)
  - Sort includes according to

     1. own header
     2. headers from the same library
     3. specific third party library (e.g., hdf5, and eigen)
     4. general third party/standard library


Code Documentation:

  - Use doxygen docstrings. In particular
  
    - Write an explicit `@brief` for method documentations. Continue the detailed description in the next line.
    - Always start a line with a capital letter and end with a dot.
    - The plural of classes, objects etc. should be denoted with a `%` sign between class name and plural s, e.g., `Household%s`. This is in order to visualize it correctly and provide a link on the doxygen page.
    - Use `[in]`, `[out]`, or `[in, out]` after `@param` in order to clarify if parameters are used as input, output or in- and output parameters.
    - To reference to enums put a # sign before the name.
    - Please also provide a description for member variables; use ``///< DESCRIPTION`  or `/**< DESCRIPTION */`` for two lines. Keep it short.


Mandatory C++ Style Guidelines
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The style guidelines are adopted from `TiGL <https://github.com/DLR-SC/tigl>`_.


**Tabs and Indentation**

 - Use 4 spaces indentation. **Don't use tabs!**
 - Exceptions:

   - public/protected/private keywords in class definitions
   - namespaces

.. code:: cpp

    namespace mio
    {
    namespace foo
    {
    namespace bar 
    {
        /*some code*/
    } // namespace bar
    } // namespace foo
    } // namespace mio

**Definitions and Declarations**

 - Braces in new lines:

.. code:: cpp

    class Secir
    {

    private:
        double m_member;
    };

If you use several lines for a functions definition/declaration, align the function arguments horizontally:

.. code:: cpp

    ReturnCode compute_something(Arg1 arg1,
                                 Arg2 arg2,
                                 Arg3 arg3)


**Loops, If and Switch Statements**

 - space before and after condition
 - Braces in the same line

.. code:: cpp

    if (psi.size()<=2) {
        psi.clear();
    }
    else {
        double psimax = psi[psi.size()-1];
    }


    for (size_t i = 0; i < psi.size(); i++) {
        some code
    }


    switch (GetSymmetryAxis()) {
    case TIGL_X_Y_PLANE:
        return zmax - zmin;
    case TIGL_X_Z_PLANE:
        return ymax - ymin;
    }


**Automatic code formatting with clang-format**

The Clang-Format Tool can also be used to reformat the code to our style. Here are the settings that should comply to our style.

.. code::

    BasedOnStyle: LLVM
    IndentWidth: 4
    SortIncludes:    false
    ColumnLimit:     120
    AlignTrailingComments: false
    AccessModifierOffset: -4
    AlignConsecutiveAssignments: true
    ReflowComments:  false
    BraceWrapping:   
    AfterClass:    true
    AfterFunction: true
    BeforeElse: true
    BeforeCatch: true
    AfterNamespace:  true
    AfterEnum: true
    BreakBeforeBraces: "Custom"
    PointerAlignment: Left
    AllowShortFunctionsOnASingleLine: false
    NamespaceIndentation: None
    BreakConstructorInitializersBeforeComma: true
    AlwaysBreakTemplateDeclarations: Yes
    AllowShortLambdasOnASingleLine: Empty


These settings are set in the file ``.clang-format`` in the root directory of the repository. 

**Using clang-format with either Qt, Visual Studio Code, or VSCodium**

The Beautifier plugin shipped with QtCreator supports clang-format (help could also be provided by https://www.vikingsoftware.com/using-clang-format-with-qtcreator/ ), so you will be able to automatically format your code. For Visual Studio Code, install the Clang-format extension and add the lines:

.. code:: 

    "editor.formatOnSave": true,
    "clang-format.executable": "...path...to...clang-format-executable",

to your settings.json and store the above code formatting rules in a file named ``.clang-format`` in the working directory of VSCode.

Note: The clang-format provided by default in Debian/Ubuntu is quite old and with our style file the issue

.. code:: bash

    YAML:21:34: error: invalid boolean
    AlwaysBreakTemplateDeclarations: Yes
                                    ^~~
    Error reading PATH/.clang-format: Invalid argument


might appear. In that case, update ``clang-format`` or install a newer version (e.g. ``clang-format-10``) manually and point to its executable.


Python coding guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please follow the `PEP 8 -- Style Guide for Python. <https://www.python.org/dev/peps/pep-0008/>`_


**Note on maximum line length**

If using autopep8, e.g., of the Python plugin for Visual Studio Code or VSCodium, maximum length might not be correctly applied. In this case, add

.. code::

    "python.formatting.autopep8Args": ["--max-line-length", "79", "--experimental"]

to your corresponding ``settings.json``.


**Docstrings**

Docstrings in Python should be added for every function, as detailed in the C++ coding guidelines. However, the syntax is slightly different than for C++ code. An overview and examples can be found at https://sphinx-rtd-tutorial.readthedocs.io/en/latest/docstrings.html . 

Figure colors and settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to ensure that figures in the documentation and in the code have a consistent look, we use the following settings:

**Default color scheme**

- For figures in the documentation, we usually use the `matplotlib <https://matplotlib.org/>`_ library. 
- The default color cycle is set to the `Set1 <https://matplotlib.org/stable/tutorials/colors/colormaps.html#Qualitative>`_ colormap.

**Colorblind-friendly alternatives**

For better accessibility and when creating figures with many categories, consider using colorblind-friendly alternatives:

- Use the `tab10 <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_ colormap for up to 10 distinct categories
- For sequential data, prefer `viridis <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_, `plasma`, or `cividis` colormaps
- For diverging data, use `RdBu <https://matplotlib.org/stable/tutorials/colors/colormaps.html>`_ or `RdYlBu` colormaps
- Avoid using red-green color combinations without additional visual cues (patterns, shapes, etc.)

**General figure guidelines**

- Use consistent font sizes across all figures (typically 10-12pt for labels, 8-10pt for tick labels)
- Ensure sufficient contrast between colors and background
- Add appropriate legends and axis labels with units
- For line plots with multiple series, vary both color and line style (solid, dashed, dotted) for better distinction
- When possible, test figures with a colorblind simulator to ensure accessibility

Git workflow
----------------------

General
~~~~~~~~~~~~

- There is a main but no release or develop branch. The main branch is always stable and the latest release can be found as a tagged commit on the main branch. Stable means that all tests pass.
- All actual work is done in task branches, regardless of whether it's a feature, a bugfix, or a performance analysis.
- Task branches are generally created from the main branch.
- Please **never rebase** your branches, **always use merging** (with the main or other changes) such that committed changes can be followed in the history. There will be a squashed commit when the changes are added to the main branch.
- The name of a task branch satisfies the following template: ``issueId-issueName``.
- Each commit must have a meaningful commit message.
- In general, we should try to keep the branches working. However, if you need to commit a non-working change, please begin the commit message with ``[ci skip] non-working`` for the following reasons:
  
  - Nobody attempts to checkout this commit with the assumption it would work.
  - ``[ci skip]`` prevents the CI from running.

- If we release a new version of the software, we create a tag for the version on the main branch.
- Please keep all issue-related communication within the issue or pull request.

Software Development in Sprints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The software development process is inspired by `Scrum <https://en.wikipedia.org/wiki/Scrum_(software_development)>`_ and the development of the core developers is organized in sprints. The rules below only partially apply to external (non-core) contributors.

**General**


- A sprint is a temporally limited cycle of a fixed time, in our case **three** weeks.
- The scope of work will be defined in a sprint meeting where work is related to issues.
- MEmilio-related issues are categorized in three different classes: agent-based modeling, equation-based modeling and MEmilio: data, tools and more. If a clear categorization is not possible, issues may be assigned to more than one class.
- Sprints are organized via the new GitHub Project boards: https://github.com/DLR-SC/memilio/projects

**Procedure**

- At the latest in the morning *before* every sprint meeting, all developers are encouraged to think about which issues should be processed in the upcoming sprint, regardless of whether those issues are tasks for oneself or someone else. Therefore, those issues are marked with the upcoming project. Every developer should put the issues they want to work on or request others to work on in the "SprintBacklog" and attribute them with the next sprint number.
- Shortly before the meeting, **every** developer should already look at the project issues and think about the time needed for realization.
- In the meeting, we go through the different issues and clarify questions and comments.

**New Tasks**

- For every single programming task, bug report, discussion item etc., open a new issue.
- Every issue should contain a detailed description of the task and subtasks understandable for all developers in the project.
- A new issue has no status label. Additional labels should be appended, see the label list below. At this point, it is not necessary to assign it to someone.
- Every issue should be tagged with at least one of the projects, if possible.
- Tasks (issues) which are attributed to a sprint are tracked in an `issue board <https://github.com/DLR-SC/memilio/projects>`_ found under "Projects".

**Working on an Issue**


- When you start working on an issue, make sure it is attributed to the current sprint.
- Then, assign it to yourself and move it into the column "In Progress" or change the label to ``status::in progress``. If code changes are involved, create a branch. If you are working with a partner, mention the partner in the issue description. The assignee is responsible for the ticket.
- You should only work on one ticket at a time. It shouldn't happen that two tickets are "In Progress" at the same time.
- If you completed the issue, set the pull request to "Ready for Review". Check that all coding requirements of the author (automatically added as checkboxes) are met. Assign the pull request to the person who should review your work and move the issue into the column "in review" or change the status to ``status::in review``.

**Review**


- The task of the reviewer is to check the code for correctness, compliance with the coding guidelines, completeness, quality, etc.
- If the review has finished and there are complaints, the issue is moved to ``status::in progress`` again, reassigned to the original assignee, and the threads must be resolved. Add the ``WIP`` tag to the merge request again.
- If the reviewer approves the work, the new code can be merged and the issue is automatically "Closed".
- The reviewer is allowed to make small changes in the comments or documentation, e.g., remove typos. Also, small changes such as adding/deleting spaces or empty lines can be made directly by the reviewer to speed up the process.
- For all other changes concerning the code and its functionality, the reviewer has to write a comment or make a suggestion which then goes back to the developer (see above).

**Authors and Contributions**


To honor original authors as well as reviewers and their suggestions, reviewers should be added as co-authors when merging a pull request. To do so, add the following line(s) at the end of the commit message:

.. code-block:: text

    COMMIT_MSG

    Co-authored-by: NAME <ADDRESS@XYZ.COM>

**Label List**


The full list of labels that should be used to identify issues can be found at: https://github.com/DLR-SC/memilio/labels


Documentation
--------------------

The documentation uses `Sphinx <https://www.sphinx-doc.org/en/master/>`_ and is written in reStructuredText, that uses a 
slightly different syntax than Markdown. A documentation can be found `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_.
This online documentation is generated using `ReadTheDocs <https://readthedocs.org/>`_ and is automatically updated when 
a pull request is merged into the main branch. Thus, we recommend building the documentation locally to test changes.


Please make sure to have a working python environment with a python version that is compatible with 
our :doc:`memilio-python packages <python/python_packages>` as well as 
all packages listed in ``docs/requirements.txt`` and `doxygen <https://doxygen.nl/>`_ installed.

First generate the doxygen output by running 

.. code-block:: bash

    cd docs
    doxygen


In the ``docs/Doxyfile`` (line 736), you can change for which folders the doxygen output should be generated. For faster 
build times while testing we recommend to only use e.g. ``../cpp/models/abm``. PLEASE don't commit this change!

Then sphinx can be used to build the documentation:

.. code-block:: bash

    cd docs
    make html # sphinx-build source html

The generated documentation can be found in ``docs/build/html`` (``docs/source/html`` if built without make).

For the documentation, please keep in mind that it is written in reStructuredText (RST) and uses a slightly different syntax than Markdown. A documentation can be found at `<https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html>`_.
