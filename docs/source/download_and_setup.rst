Download and setup
==================

MEmilio consists of a core library written in C++, including several models and examples, as well as several Python
packages

We have structured this documentation to guide you step-by-step through the installation and usage process.
If you still need help, feel free to `contact us <mailto:Martin.Kuehn@DLR.de>`__ or open an issue at `GitHub <https://github.com/SciCompMod/memilio/issues>`__ 
and highlight @mknaranja and @HenrZu such that we can assist you as best as we can.

There are two main ways to set up MEmilio on your computer or on a remote cluster or supercomputer,
depending on what you want to do:

1.  **Using the Python packages:** This is the recommended path for many users not familiar with C++.
    With the :doc:`memilio-simulation <m-simulation>` package you can run simulations through Python bindings.
   
    If you have a working Python installation, head to :ref:`python-package-installation` and install the package from
    PyPI (Option 1).
    For setting up Python, or for installing our other Python packages, continue here for the setup instructions.

2.  **Directly building the C++ Core:** This is for developers who want to modify the functionality,
    contribute new models etc. by running C++ code directly.

In addition, we provide several Python packages to download epidemiological data or create plots from Python.

Below, we will give you a step-by-step guide for both methods. If you are new to MEmilio and more familiar with Python, 
Julia, or R than with C++, we recommend starting with the Python packages, as they provide an easy access to simulate 
infection dynamics models and gain experience with MEmilio.

.. _general-setup-section:

Setup: Required tools
---------------------

Before you can install MEmilio, you need to install some common development tools. 

*   **Python:** Required for the Python packages.

    MEmilio is tested daily with `Python <https://www.python.org/>`__ 3.9 and 3.13.
    While other versions may also work, we recommend using the latest release of either of these. 
    
    .. dropdown:: :fa:`gears` Python installation help

        You can check whether Python is installed by opening a terminal and running

        .. code-block:: console

            python --version

        or 

        .. code-block:: console

            python3 --version

        to see which version is installed. If no installation can be found it will say "command not found" or something similar.
        You can also type ``python`` and hit the tab key once or twice, it may auto-complete to a Python command 
        including the minor version, e.g. ``python3.13``. Add the ``--version`` flag and hit enter.

        If one of these commands reports a supported version, you are done here!
        Do make note of which command works for your setup, as we will simply use ``python`` in this documentation.

        Otherwise, you need to install python first. Here are some platform specific hints:

        *   On **Linux**, all common distributions come preinstalled with a version of python 3.
            It should work out of the box.

            Optionally, you can use tools like `pyenv <https://github.com/pyenv/pyenv>`__,
            `mise <https://mise.jdx.dev/>`__, and `asdf <https://asdf-vm.com/>`__ to install an additional python 
            version.

            Be careful when installing a different version with your distribution's package manager,
            as replacing the main python installation may break your system.

        *   On **MacOS**, you can use `homebrew <https://brew.sh/>`__ to install a specific python version,
            or use the official installer from https://www.python.org/downloads/.
            A guide can be found `here <https://docs.python.org/3/using/mac.html>`__.

        *   On **Windows**, you can use `winget <https://learn.microsoft.com/en-us/windows/package-manager/winget/>`__
            to install a specific python version, or use the Python Install Manager from the
            `Microsoft Store <https://apps.microsoft.com/detail/9nq7512cxl7t>`__.
            A guide can be found `here <https://docs.python.org/3/using/windows.html>`__.

        *   On a cluster or supercomputer, python may be included as a module (e.g. through spack) that has to be
            loaded first. Check the official documentation of that system for details.

        If you are using another python management solution like conda, check out their documentation on how to set up
        an environment.

    .. note::

        In this documentation, we will simply use ``python`` as a command, but depending on your installation method
        and platform you may need to substitute this by, e.g., ``python3``, ``python3.11``, ``.\python3.13.exe``, or
        similar.

*   **C++ Compiler and CMake:**

    *   **Windows:** The easiest way is to install **Visual Studio Community**.
        This includes a C++ compiler, CMake, and Git all in one.
    
    *   **macOS:** One option is installing the **Xcode Command Line Tools** by running
        ``xcode-select --install`` in your terminal.
    
    *   **Linux:** On Linux, essential build tools and CMake might be preinstalled. Otherwise, on Debian/Ubuntu,
        you can install them by running ``sudo apt-get install cmake gcc g++`` in your terminal.

*   **Git:** This is a version control system used to download the project's source code.

    On Linux and MacOS, Git is often preinstalled, but not on Windows.
    You can check by opening a terminal and running ``git --version``.
    Instructions on how to download and install it can be found on `git-scm.com <https://git-scm.com/install>`__.

.. _general-download-section:

Download the MEmilio source code
--------------------------------

For C++ developement and for building python from source, it is required to download MEmilio.

You can get the newest version (or a specific release) of MEmilio from our
`GitHub repository <https://github.com/SciCompMod/memilio>`__, preferably by using Git:

Open a terminal and navigate to a directory of your choice, then run

.. code:: bash

    git clone https://github.com/SciCompMod/memilio

.. dropdown:: :fa:`gears` A Note on HTTPS vs. SSH

    The ``git clone`` command above uses an **HTTPS** URL. This is the simplest method and works perfectly for
    downloading the code.

    However, if you plan to contribute code back to the project (i.e., "push" your changes), we recommend using **SSH**.
    This requires a GitHub Account as well as proper authentication. To set this up, you can follow
    `GitHub's official guide on adding an SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__.

    Then, use the URL ``git@github.com:SciCompMod/memilio.git`` for cloning instead. You can also change to the ssh
    address later using

    .. code:: bash

        git remote set-url origin git@github.com:SciCompMod/memilio.git

This will create a new directory called ``memilio`` on your computer, containing the MEmilio project.
In your terminal, change into this directory:

.. code-block:: console

   cd memilio

Next steps
----------

From here, you can

* :doc:`continue with the C++ build instructions <cpp/installation>` or

* :ref:`start building the Python packages from source (Option 2) <python-package-installation>`.

Further questions
-----------------
If you have any further questions, please take a look at our :doc:`faq` and feel free to contact us via
`e-mail <mailto:Martin.Kuehn@DLR.de>`__ or open an issue or discussion on
`GitHub <https://github.com/SciCompMod/memilio/discussions>`__.
