Getting Started
===============

Overview
-------------

.. note:: This project is under active development.


MEmilio is an extensive framework for tasks around infectious disease modeling. It supports a multitude of :ref:`model <model-faq>` types 
including :doc:`equation-based<cpp/aggregated_models>`, :doc:`agent-based <cpp/individual_models>`, 
and :doc:`hybrid graph-ODE-based models <cpp/graph_metapop>`. It furthermore provides ready-to-use tools for data integration and visualizations. 
Among the equation-based models, we provide models based on :doc:`ordinary differential equations <cpp/ode>`,
:doc:`the linear chain trick (LCT), <cpp/lct>` and a recent :doc:`generalized LCT <cpp/glct>`, :doc:`integro-differential equations <cpp/ide>` 
and :doc:`stochastic differential equations <cpp/sde>`. With simple definitions, models can be spatially or demographically resolved.

The MEmilio framework is written in two languages: C++ and Python. 

- The C++ backend contains efficient and optimized model implementations that further use parallelization to speed up execution and reduce waiting times.
- Python is used for data acquisition, plotting, and machine-learning models.
- We, furthermore, provide Python interfaces to selected models (implemented in C++) to allow the use and study of advanced models by users less experience in programming or computer science.

For more details on using models implemented in C++ directly, see the sections on :doc:`model usage <cpp/model_usage>`.
For more details on implementing new infection dynamics models that could then be combined with, e.g., our mobility patterns, see :doc:`model creation <cpp/model_creation>`.

If you prefer using Python to call or run our models, you can use our :doc:`memilio-simulation <python/memilio_simulation>` package to run simulations.
The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download and structure important data such 
as infection or mobility data. More about this and our other Python packages can be found in the :doc:`Python Interface Section <python/python_packages>` 
of this documentation.

A few things are not represented in this documentation, but part of the `github repository <https://github.com/SciCompMod/memilio>`_. 
In the `data <https://github.com/SciCompMod/memilio/tree/main/data>`_ folder you can find some regularly used data 
for simulations of a pathogen's spread, currently mostly for Germany. 


Usage
-----------------

.. _installation:

Installation
~~~~~~~~~~~~

There are two main ways to set up MEmilio on your computer or on a remote cluster or supercomputer, depending on what you want to do:

1. **Using the Python packages:** This is the recommended path for many users not familiar with C++. Here, you can run simulations using python bindings.
2. **Directly building the C++ Core:** This is for developers who want to modify the functionality, contribute new models etc. by running C++ code directly.

In addition, we provide several Python packages to download epidemiological data or create plots from Python.

Below, we will give you a step-by-step guide for both methods. If you are new to MEmilio and more familiar with Python, Julia, or R than with C++, we recommend starting with the Python packages, 
as they provide an easy access to simulate infection dynamics models from and collect experiences with MEmilio.

Required Tools
*****************

Before you can install MEmilio, you need to install some common development tools. 

*   **Git:** This is a version control system used to download the project's source code.

    *   **Windows:** By default, Git is not installed. Download and install it from `git-scm.com <https://git-scm.com/downloads/win>`_.
    *   **macOS & Linux:** Git is usually preinstalled. You can check by opening a terminal and typing ``git --version``.

*   **Python:** Required for the Python packages.

    *   MEmilio is tested daily with Python 3.8 and 3.11. While other versions might also work, we recommend installing the latest version tested daily from the official website `python.org <https://www.python.org/>`_.

*   **C++ Compiler and CMake:**

    *   **Windows:** The easiest way is to install **Visual Studio Community**. This includes a C++ compiler, CMake, and Git all in one.
    *   **macOS:** One option is installing the **Xcode Command Line Tools** by running ``xcode-select --install`` in your terminal.
    *   **Linux:** On Linux, essential build tools and CMake might be preinstalled. Otherwise, on Debian/Ubuntu, you could execute the installation by running ``sudo apt-get install cmake gcc g++`` in your terminal.

Step 1: Download the MEmilio Source Code
*****************************************

Once the required tools are installed, open a terminal and download the MEmilio code with this command:

.. code-block:: console

   git clone https://github.com/SciCompMod/memilio.git

This command copies the entire MEmilio project into a new folder named ``memilio`` on your computer. 

.. note:: A Quick Note on HTTPS vs. SSH

   The ``git clone`` command above uses an **HTTPS** URL. This is the simplest method and works perfectly for downloading the code.

   However, if you plan to contribute code back to the project (i.e., "push" your changes), we recommend using **SSH**. To set this up, you can follow `GitHub's official guide on adding an SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`_.


Now, navigate into that folder:

.. code-block:: console

   cd memilio

From here, choose one of the following options.

Option A: Installing the Python Packages (Recommended for nonexperienced users or for data download and visualizations)
****************************************************

You can run simulations, download data, or create plots, by only installing our Python packages.

1.  Navigate to the directory containing our Python code:

    .. code-block:: console

       cd pycode

2.  To install the simulation package ``memilio-simulation``, from here you can do:

    .. code-block:: console

       cd memilio-simulation
       pip install -e .

3.  For afterwards installing the ``memilio-epidata`` package for data downloading and handling, run:

    .. code-block:: console

       cd ..  # Go back to the pycode directory
       cd memilio-epidata
       pip install -e .

.. tip:: For Contributors: Installing development packages

   The ``-e`` flag installs the package in a mode, which links the installation to your local source code folder.

   If you plan to contribute to MEmilio, you can also install all the necessary development dependencies by adding ``[dev]`` to the command:

   .. code-block:: console

      pip install -e .[dev]

   For regular use, the simple ``pip install -e .`` is sufficient.

To install other packages, see the items below *Python Interface* in the menu on the left hand side.

Option B: Building the C++ Core (Advanced)
****************************************

For experienced developers and C++ programmers, we offer the C++ backend to fully benefit from all functionality and parallel performance.

1.  Navigate to the C++ source code directory:

    .. code-block:: console

       cd cpp

2.  Create a separate directory for the build files.

    .. code-block:: console

       mkdir build && cd build

3.  Run CMake. This tool prepares the project for compilation on your specific system.

    .. code-block:: console

       cmake ..

4.  Compile the code and create the executables.

    .. code-block:: console

       cmake --build .

For more detailed instructions, help with errors, and a list of compile options, please see the full :doc:`C++ Installation Guide <cpp/installation>`.

Running simulations
~~~~~~~~~~~~~~~~~~~~~
You can run simulations either via the C++ interface where they are originally implemented or via the python bindings. 
For the C++ Interface, you can find explanations of the models as well as guides on their usage in the :doc:`C++ model usage <cpp/model_usage>` section.
In short, the executables for different model instantiations are build as described above and can be run via 

.. code-block:: console

   ./cpp/build/bin/<example_name>


Out of the box this works for all examples in the ``cpp/examples`` folder of our `github repository <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`_,
that do not depend on user-provided external libraries. 
Additional explanations for our models are linked at the corresponding sites of this documentation.

Simulations used in publications
********************************
For simulations used in publications, we maintain a separate repository: 
`memilio-simulations <https://github.com/SciCompMod/memilio-simulations>`_. 
This repository contains simulations organized in separate folders, each with the specific version of MEmilio 
used for the published results. This ensures that simulation results can be easily reproduced.

The repository also includes additional scripts for plotting, data gathering, and pre-/post-processing 
that were used in publications.


Loading data
~~~~~~~~~~~~~~~~~~~~~
The :doc:`memilio-epidata <python/memilio_epidata>` package provides tools to download epidemiological relevant datasets. Some 
datasets like contact matrices for Germany are also included in the ``data`` folder of the `github repository <https://github.com/SciCompMod/memilio/tree/main/data>`_ and 
school holidays (for Germany) are directly included in the `C++ code <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/geography/holiday_data.ipp>`_.  


Creating new models
~~~~~~~~~~~~~~~~~~~~~

If you want to create new models, you can do so via the C++ interface. For this, we recommend to have a look at 
the :doc:`C++ model creation <cpp/model_creation>` section of this documentation.


Visualizations
~~~~~~~~~~~~~~~~~~~~~

For visualizations, we provide our :doc:`python package MEmilio-plot <python/memilio_plot>`. Apart from that, we have 
collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`_. 
For the latter, no regular testing is conducted. If you encounter errors, please `contact us <mailto:Martin.Kuehn@DLR.de>`.

Further questions
~~~~~~~~~~~~~~~~~~~~~
If you have any further questions, please take a look at our :doc:`faq` and feel free to contact us via `e-mail <mailto:Martin.Kuehn@DLR.de>` or open an issue or discusion on `github <https://github.com/SciCompMod/memilio/>`
