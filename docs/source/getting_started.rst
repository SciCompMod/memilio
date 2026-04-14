Getting started
===============

Overview
--------

.. note:: This project is under active development with approximately 10 full time developers. Any content provided in the repository has been extensively reviewed and can be considered stable. Regular function extensions can be expected. New features are, if possible, integrated in a backward-stable manner but interfaces might change over time.


MEmilio is an extensive framework for tasks around infectious disease modeling. It supports a multitude of :ref:`model <model-faq>` types 
including :doc:`equation-based or compartmental<cpp/aggregated_models>`, :doc:`agent-based <cpp/individual_models>`, 
and :doc:`hybrid graph-ODE-based models <cpp/graph_metapop>`. It furthermore provides ready-to-use tools for data integration and visualizations. 
Among the equation-based models, we provide models based on :doc:`ordinary differential equations <cpp/ode>`,
:doc:`the linear chain trick (LCT), <cpp/lct>` and a recent :doc:`generalized LCT <cpp/glct>`, :doc:`integro-differential equations <cpp/ide>` 
and :doc:`stochastic differential equations <cpp/sde>`. With simple definitions, models can be spatially or demographically resolved.

The MEmilio framework is written in two languages: C++ and Python. 

- The C++ backend contains efficient and optimized model implementations that further use parallelization to speed up execution and reduce waiting times.
- Python is used for data acquisition, plotting, and machine-learning models.
- We, furthermore, provide Python interfaces to selected models (implemented in C++) to allow the use and study of advanced models by users with less experience in programming or computer science.

For more details on using models implemented in C++ directly, see the sections on :doc:`model usage <cpp/model_usage>`.
For more details on implementing new infection dynamics models that could then be combined with, e.g., our mobility patterns, see :doc:`model creation <cpp/model_creation>`.

If you prefer using Python to call or run our models, you can use our :doc:`memilio-simulation <python/m-simulation>` package to run simulations.
The :doc:`memilio-epidata <python/m-epidata>` package provides tools to download and structure important data such 
as infection or mobility data. More about this and our other Python packages can be found in the :doc:`Python Interface Section <python/python_packages>` 
of this documentation.

A few things are not represented in this documentation, but are part of the `GitHub repository <https://github.com/SciCompMod/memilio>`__. 
In the `data <https://github.com/SciCompMod/memilio/tree/main/data>`_ folder you can find some regularly used data 
for simulations of a pathogen's spread, currently mostly for Germany. 

Why to use MEmilio
------------------

In computational epidemiology and infectious disease dynamics, models are often implemented in Python or R. However, this approach often limits the possibility to build large-scale models including an advanced level of detail, e.g., in demography, spatial resolution, or even individual immunity or to run many simulations in a short time frame. 
MEmilio addresses this challenge by providing a high-performance framework implemented in C++ that allows for large-scale modeling in short time frames to be used in research, policy advice, and education.

In the following figure, we representatively show an excerpt of Fig. 6 of `Bicker et al. (2026), DOI: 10.48550/arXiv.2602.11381 <https://doi.org/10.48550/arXiv.2602.11381>`_ showing the performance of population and metapopulation models implemented in R and in C++ in MEmilio. While for large numbers of regions, the R model based on the C-implemented routine desolve comes close to MEmilio's C++ routine performance, both interfaces (C++ and Python) of MEmilio realize significant speedups (factor 100 and more) for most applications.

.. image:: http://martinkuehn.eu/research/images/speedup_memilio.png
   :alt: Performance of population and metapopulation models implemented in R and in C++ in MEmilio.
   :width: 100%

The use of a particular model is generally driven by the research question at hand. The distinction of MEmilio is the provision of a wide range of models, from simple compartmental models to complex integro-differential and agent-based models, allowing users to select the most appropriate model for their specific needs.

:doc:`Aggregated models<cpp/aggregated_models>` are suitable for scenarios where population-level dynamics are of interest. They are computationally efficient and can be used for quick assessments or when data is limited. In our implementations, these models can easily be extended to address research questions that involve demographic dimensions such as age.

Standard models based on :doc:`ordinary differential equations (ODE) <cpp/ode>` allow the simplest description of population-level infection dynamics. However, these also implicitly assume exponentially distributed stay times. If the data suggest these to be unrealistic, ODE-based models with the :doc:`linear chain trick (LCT) <cpp/lct>` can be used. With the LCT, Gamma, or more precisely, Erlang distributions can be adapted to expected stay time and variance. We also offer a first implementation of a :doc:`generalized LCT <cpp/glct>` where more distributions can be approximated. For full flexibility, :doc:`integro-differential equation-based (IDE) models <cpp/ide>` can be used. These allow for arbitrary stay time distributions and are thus suitable when a more realistic timing of disease progression is crucial for the research question. However, IDE models are computationally more expensive than ODE-based models.

For spatio-temporal dynamics, :doc:`graph-based metapopulation models <cpp/graph_metapop>`, which leverage :doc:`ODE-based models<cpp/ode>`, are a good compromise between level of detail and computational effort. They allow for the incorporation of mobility patterns and spatial heterogeneity, making them suitable for studying the spread of diseases across different regions. They can also be used to, e.g., consider different intervention strategies in different regions or to study the effect of mobility restrictions. To study the effect of local interventions, we provide the option to implement pre-defined and dynamic NPIs which automatically enforce interventions on a local and regional level when an incidence threshold criteria by the user is exceeded. For details, see, e.g.,

*   :doc:`ODE-based SECIR model <cpp/models/osecir>` for early-phase epidemics or full immunity cases
*   :doc:`ODE-based SECIRVVS model <cpp/models/osecirvvs>` for early- or mid-phase epidemics with three immunity layers
*   :doc:`ODE-based SECIRTS model <cpp/models/osecirts>` for mid- or late-phase epidemics with three immunity layers and waning immunity

For Python, please see, e.g., :doc:`ODE-based SECIRTS model <python/m-simulation_model_usage>`.

When individual-level interactions and heterogeneity are crucial, :doc:`individual-based models <cpp/individual_models>` provide a detailed representation of disease dynamics. These models can capture complex behaviors and interactions, making them valuable for understanding transmission dynamics in specific settings. Individual-based models are computationally intensive but offer unparalleled detail for certain research questions such as in-household transmission or vaccination and testing strategies targeting individuals that satisfy specific properties with respect to age, previous infections, immunity levels, or particular workplaces. The most versatile individual-based model in MEmilio is the :doc:`(mobility-based) agent-based model <cpp/mobility_based_abm>`.

A quick tour through MEmilio
-----------------------------

While MEmilio harmonizes much of its structures across all model classes, we first distinguish between :doc:`compartmental or aggregated models<cpp/aggregated_models>` based on ODEs (ordinary differential equations) without and with Linear Chain Trick, IDEs (integro-differential equations), and SDEs (stochastic differential equations) and :doc:`Agent-based models<cpp/individual_models>`. The following subsections give a brief overview on essential functionality, each presented in a particular and function-specific tutorial. While several tutorials build on previous tutorials, experienced users or users interested in other models might also want jump to later parts of this walkthrough guide. 

An additional overview on MEmilio's elementary model structure is given by the following figure.

.. image:: http://martinkuehn.eu/research/images/model_structure.png
   :alt: Overview on MEmilio's model structure.
   :width: 100%

MEmilio benefits from a harmonized description of its models in infection states and parameters, and, potentially, a list of flows between the compartments; see the following figure for a motivation. All models derive their infection states from a flexible and simple list of InfectionStates. For FlowModels (see below for an explanation), particular transitions are defined evenly flexible as a list of flows between the states. Parameters are also generally defined in an identical fashion. 

.. image:: http://martinkuehn.eu/research/images/uniform.png
   :alt: MEmilio's uniform model description.
   :width: 100%

Below we guide you through several tutorials on using MEmilio's models through its Python interface. More experience users might directly start with the `Python exercises <https://github.com/SciCompMod/memilio-tutorials/tree/main/exercises>`_ which are derived versions from the tutorials or with `tutorial and exercises in C++ <https://github.com/SciCompMod/memilio-tutorials/tree/main/cpp-tutorials>`_. For more advanced aggregated models using the Linear Chain Trick or IDE-formulations, we currently only provide tutorials and exercises in C++. For the individual- or agent-based model (ABM), we currently only provide `ABM tutorials in C++ <https://github.com/SciCompMod/memilio-tutorials/tree/main/cpp-tutorials/abm>`_. 


Simple compartmental models
****************************

Most of MEmilio's compartmental or aggregated models share the same interface derived from a high-level **CompartmentalModel** (see above). It defines the fundamental structure for epidemiological models with compartments (e.g., SEIR, SECIR, SIRS, etc.).

In `Tutorial 01 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial01.py>`_, we show how to set up and simulate a simple setting for our :doc:`ODE-SECIR model <cpp/models/osecir>`. The result of the tutorial is a figure of a well-known epidemic outcome.

.. image:: http://martinkuehn.eu/research/images/tutorial01.png
   :alt: A well-known epidemic outcome as a result of Tutorial 01.
   :width: 100%


Flows between compartments
**************************

Often, modelers might be interested not only in the estimated number of individuals in a state of the disease but also in the number of recent or current transitions between different states such as the number of new hospitalizations. As modelers could introduce additional compartments only following those transitions or do complex post-processing, MEmilio directly computes all transitions between compartments by default. This is realized through MEmilio's **FlowModel** structure which is a still generic but refined specification of the **CompartmentalModel**. Through an optimized backend, the overhead for computing transitions (i.e. flows) and compartmental values is less than 10 %.

In `Tutorial 02 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial02.py>`_, we show how to obtain the numbers of newly symptomatic and hospitalized individuals for our :doc:`ODE-SECIR model <cpp/models/osecir>`. The result of the tutorial is shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/tutorial02.png
   :alt: Newly symptomatic and hospitalized individuals as obtained from Tutorial 02.
   :width: 100%


Demography and contact structures
*********************************

As motivated in the following figure, MEmilio's models are implemented in a way that they can be stratified by age groups (or other dimensions such as sex) in a single line.

.. image:: http://martinkuehn.eu/research/images/contacts.png
   :alt: Module for flexible demographic stratification by age groups.
   :width: 100%

In `Tutorial 05 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial05.py>`_, we show how to distinguish individuals of three different age groups by their susceptibility with respect to severe and critical infections and simulate outcomes for our :doc:`ODE-SECIR model <cpp/models/osecir>`. The result of the tutorial is shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/tutorial05.png
   :alt: Different epidemic curves for six different age groups.
   :width: 100%

Metapopulation and mobility
***************************

As motivated in the following figure, MEmilio's aggregated models can be extended to metapopulation models by using a graph structure.

.. image:: http://martinkuehn.eu/research/images/mobility.png
   :alt: Module for flexible spatial resolution in metapopulation models.
   :width: 100%

In `Tutorial 07 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial07.py>`_, we show how an epidemic with our :doc:`ODE-SECIR model <cpp/models/osecir>` evolves with a delay between two different spatial entities. The result of the tutorial is shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/tutorial07.png
   :alt: Delayed epidemic spreading through metapopulation coupling of two regions.
   :width: 100%

Fixed time-point interventions
******************************

In order to control and mitigate epidemic developments, MEmilio provides the ability to introduce non-pharmaceutical interventions (NPIs) or measures as `Dampings` to the contact frequencies. In `Tutorial 03 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial03.py>`_, we show how an epidemic with our :doc:`ODE-SECIR model <cpp/models/osecir>` can be first mitigated before a reopening event takes place. The result of the tutorial is shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/tutorial03.png
   :alt: Changed epidemic outcome through interventions at fixed time points.
   :width: 100%


Location-specific interventions
*******************************

Often interventions are targeted to specific types of locations such as schools, workplaces, or social gatherings. In order to most realistically model contact structures and NPIs across different locations, MEmilio uses simple and flexible lists of contact locations. In `Tutorial 10 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial10.py>`_, we explain with our :doc:`ODE-SECIR model <cpp/models/osecir>` how contact structures can be stratified by locations and NPIs implemented in a location-specific way. The result of the tutorial is shown in the following figure.

.. image:: http://martinkuehn.eu/research/images/tutorial10.png
   :alt: Changed epidemic outcome through interventions at specific locations.
   :width: 100%

Dynamic interventions
*********************

Eventually, NPIs might often be bound to a threshold or criterion upon which its get activated, e.g., the number of new symptomatic (here, reported) infections over the last days. In order to allow dynamic, threshold-dependent NPIs, MEmilio implements a structure denoted `DynamicNPIs`. In `Tutorial 11 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial11.py>`_, we explain with our :doc:`ODE-SECIR model <cpp/models/osecir>` how to set up and simulate dynamic interventions based on symptomatic infections. The result of the tutorial is shown in the following figure. 

Note that the DynamicNPI feature is currently fixed to interventions based on symptomatic infections but if you are interested in using it for other applications, please get in touch with us, as the change could be done by us in very short time.

.. image:: http://martinkuehn.eu/research/images/tutorial11.png
   :alt: Changed epidemic outcome through dynamically activated interventions.
   :width: 100%

Fitting MEmilio's models
************************

As parameter inference is a research topic of its own, MEmilio does not provide methods for parameter inference but instead provides well designed interfaces to established tools and packages dedicated to model calibration and parameter inference.

`Tutorial 04 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial04.py>`_ and `Tutorial 06 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial06.py>`_, we introduce usage of Approximate Bayesian Computation (ABC) with MEmilio and `pyABC <https://pyabc.readthedocs.io/en/latest/>`_ for likelihood-free inference. 

The result of Tutorial 04 are the projections of the calibrated model using pyABC:

.. image:: http://martinkuehn.eu/research/images/tutorial04.png
   :alt: Projections of the calibrated model using pyABC.
   :width: 100%

The result of Tutorial 06 are the posterior distributions for the model parameters using pyABC:

.. image:: http://martinkuehn.eu/research/images/tutorial06.png
   :alt: Posterior distributions for the model parameters using pyABC.
   :width: 100%

In `Tutorial 09 <https://github.com/SciCompMod/memilio-tutorials/blob/main/tutorial09.py>`_ we use `Bayesflow <https://bayesflow.org/main/index.html>`_, a state of the art python library for Bayesian inference with deep learning. The result of Tutorial 09 is the region- and age-specific calibration using BayesFlow:

.. image:: http://martinkuehn.eu/research/images/tutorial09.png
   :alt: Region- and age-specific calibration using BayesFlow.
   :width: 100%

Linear Chain Trick
*******************

As among others shown in `Plötzke et al. (2026), DOI: 10.1016/j.matcom.2025.07.045 <https://doi.org/10.1016/j.matcom.2025.07.045>`_, exponentially distributed stay times can lead to (substantially) deviating peak timings and values; see also the following figure extracted from Fig. 8 of the mentioned paper.

.. image:: http://martinkuehn.eu/research/images/lct.png
   :alt: Deviating peaks with exponential distribution assumptions (ODE) versus Linear Chain Trick with Erlang/Gamma distributions.
   :width: 100%

In the `LCT Tutorial <https://github.com/SciCompMod/memilio-tutorials/blob/main/cpp-tutorials/tutorial_lct.cpp>`_, we show a minimalistic example of an LCT model within MEmilio. Parameters, contact structures, and NPIs can basically be used as in the introductions to simple ODE models. If you are interested in using LCT models with different courses of the disease, please get in touch with us.

IDE-based models
*****************

While LCT models generalize assumptions from exponential to Gamma, IDE formulations such as presented in `Wendler et al. (2026), DOI: 10.1016/j.amc.2025.129636 <https://doi.org/10.1016/j.amc.2025.129636>`_ allow the full flexibility to use any data-driven transition distribution.

In the `IDE Tutorial <https://github.com/SciCompMod/memilio-tutorials/blob/main/cpp-tutorials/tutorial_ide.cpp>`_, we show a minimalistic example of an IDE model within MEmilio. Parameters, contact structures, and NPIs can basically be used as in the introductions to simple ODE models. If you are interested in using IDE models with different courses of the disease, please get in touch with us.


Agent- or individual-based model
********************************

As motivated in Fig. 5 of `Bicker et al. (2026), DOI: 10.48550/arXiv.2602.11381 <https://doi.org/10.48550/arXiv.2602.11381>`_ and shown here below, different types of models with their implicit or flexible assumptions can lead to all type of epidemic projections.

.. image:: http://martinkuehn.eu/research/images/model_comparisons.png
   :alt: Results of infectious disease spread with different models and model assumptions.
   :width: 100%

A major advantage of agent-based models (ABMs) is the possibility to model individuals with individual properties. In the MEmilio-ABM, populations are set up straightforward with household structures. Testing and vaccination strategies can be targeted to individuals at particular locations or of selected age groups.

- In the `tutorial on households <https://github.com/SciCompMod/memilio-tutorials/blob/main/cpp-tutorials/abm/tutorial_abm_households.cpp>`_, we show how particular populations and household structures can be set up with the MEmilio-ABM.
- In the `tutorial on testing <https://github.com/SciCompMod/memilio-tutorials/blob/main/cpp-tutorials/abm/tutorial_abm_testing.cpp>`_, we show how different testing strategies with particular testing schemes can be realized.
- In the `tutorial on vaccination <https://github.com/SciCompMod/memilio-tutorials/blob/main/cpp-tutorials/abm/tutorial_abm_vaccination.cpp>`_, we show how different vaccination strategies can be realized.


Download and installation
-------------------------

The installation and use of MEmilio might look overwhelming at first due to the many features and models included. 
We have structured this documentation to guide you step-by-step through the installation and usage process.
If you still need help, feel free to `contact us <mailto:Martin.Kuehn@DLR.de>`_ or open an issue at `GitHub <https://github.com/SciCompMod/memilio/issues>`_ 
and highlight @mknaranja and @HenrZu such that we can assist you as best as we can.

.. _installation:

Installation
~~~~~~~~~~~~

There are two main ways to set up MEmilio on your computer or on a remote cluster or supercomputer, depending on what you want to do:

1. **Using the Python packages:** This is the recommended path for many users not familiar with C++. Here, you can run simulations using python bindings.
2. **Directly building the C++ Core:** This is for developers who want to modify the functionality, contribute new models etc. by running C++ code directly.

In addition, we provide several Python packages to download epidemiological data or create plots from Python.

Below, we will give you a step-by-step guide for both methods. If you are new to MEmilio and more familiar with Python, Julia, or R than with C++, we recommend starting with the Python packages, as they provide an easy access to simulate infection dynamics models from and collect experiences with MEmilio.

Required tools
**************

Before you can install MEmilio, you need to install some common development tools. 

*   **Git:** This is a version control system used to download the project's source code.

    *   **Windows:** By default, Git is not installed. Download and install it from `git-scm.com <https://git-scm.com/downloads/win>`__.
    *   **macOS & Linux:** Git is usually preinstalled. You can check by opening a terminal and typing ``git --version``.

*   **Python:** Required for the Python packages.

    *   MEmilio is tested daily with Python 3.8 and 3.12. While other versions may also work, we recommend using the latest release of either of these. You can download it from the official website `python.org <https://www.python.org/>`__.

*   **C++ Compiler and CMake:**

    *   **Windows:** The easiest way is to install **Visual Studio Community**. This includes a C++ compiler, CMake, and Git all in one.
    *   **macOS:** One option is installing the **Xcode Command Line Tools** by running ``xcode-select --install`` in your terminal.
    *   **Linux:** On Linux, essential build tools and CMake might be preinstalled. Otherwise, on Debian/Ubuntu, you could execute the installation by running ``sudo apt-get install cmake gcc g++`` in your terminal.

Step 1: Download the MEmilio source code
****************************************

Once the required tools are installed, open a terminal and download the MEmilio code with this command:

.. code-block:: console

   git clone https://github.com/SciCompMod/memilio.git

This command copies the entire MEmilio project into a new folder named ``memilio`` on your computer. 

.. note:: A Quick Note on HTTPS vs. SSH

   The ``git clone`` command above uses an **HTTPS** URL. This is the simplest method and works perfectly for
   downloading the code.

   However, if you plan to contribute code back to the project (i.e., "push" your changes), we recommend using **SSH**.
   To set this up, you can follow
   `GitHub's official guide on adding an SSH key <https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account>`__.


Now, navigate into that folder:

.. code-block:: console

   cd memilio

From here, choose one of the following options.

Option A: Installing the Python packages (Recommended for nonexperienced users or for data download and visualizations)
***********************************************************************************************************************

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

Option B: Building the C++ core (Advanced)
******************************************

For experienced developers and C++ programmers, we offer the C++ backend to fully benefit from all functionality and parallel performance.

Please see the full :doc:`C++ Build instructions <cpp/installation>` for more details and a list of compile options.

1.  Run CMake. This tool *configures* the project for compilation on your specific system. It takes around 10 seconds, depending on your internet connection as external libraries are fetched.

    .. code-block:: console

        cmake -S cpp -B cpp/build

2. Compile the code and create the executables.
Run the build command from inside your build directory. To speed up the process,
you can use the ``-j`` flag (e.g., using 4 cores):

.. code-block:: bash

   cmake --build . -j 4

.. note::
   On a standard 4-core (2024) laptop, compilation takes approximately 6 minutes.
   Upon completion, the executables are located in the ``cpp/build/bin`` directory.

    .. code-block:: console

        cmake --build cpp/build

If you want to build a specific example, you can specify it with the ``--target`` flag:

.. code-block:: console

   cmake --build . --target <example_name>

If you experience errors, feel free to contact martin.kuehn@dlr.de or open a `discussion on GitHub <https://github.com/SciCompMod/memilio/discussions>`_!

Running simulations
~~~~~~~~~~~~~~~~~~~
You can run simulations either via the C++ interface where they are originally implemented or via the python bindings. 
For the C++ Interface, you can find explanations of the models as well as guides on their usage in the :doc:`C++ model usage <cpp/model_usage>` section.
In short, the executables for different model instantiations are built as described above and can be run via 

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
~~~~~~~~~~~~
The :doc:`memilio-epidata <python/m-epidata>` package provides tools to download epidemiological relevant datasets. Some 
datasets like contact matrices for Germany are also included in the ``data`` folder of the `github repository <https://github.com/SciCompMod/memilio/tree/main/data>`_ and 
school holidays (for Germany) are directly included in the `C++ code <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/geography/holiday_data.ipp>`_.  


Creating new models
~~~~~~~~~~~~~~~~~~~

If you want to create new models, you can do so via the C++ interface. For this, we recommend to have a look at 
the :doc:`C++ model creation <cpp/model_creation>` section of this documentation.


Visualizations
~~~~~~~~~~~~~~

For visualizations, we provide our :doc:`python package MEmilio-plot <python/m-plot>`. Apart from that, we have 
collected some scripts that we used for visualizations in the `tools folder in our github repository <https://github.com/SciCompMod/memilio/tree/main/tools>`_. 
For the latter, no regular testing is conducted. If you encounter errors, please `contact us <mailto:Martin.Kuehn@DLR.de>`_.

Further questions
~~~~~~~~~~~~~~~~~
If you have any further questions, please take a look at our :doc:`faq` and feel free to contact us via
`e-mail <mailto:Martin.Kuehn@DLR.de>`_ or open an issue or discussion on `GitHub <https://github.com/SciCompMod/memilio/>`_.
