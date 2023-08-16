# MEmilio - a high performance Modular EpideMIcs simuLatIOn software #

![memilio_logo](docs/memilio-small.png)

[![CI](https://github.com/DLR-SC/memilio/actions/workflows/main.yml/badge.svg)](https://github.com/DLR-SC/memilio/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/DLR-SC/memilio/branch/main/graph/badge.svg?token=DVQXIQJHBM)](https://codecov.io/gh/DLR-SC/memilio)

MEmilio is a common project between the Institute for Software Technology of the German Aerospace Center (DLR) and the department of Systems Immunology (SIMM) of the Helmholtz Center for Infection Research (HZI). This project will bring cutting edge and compute intensive epidemiological models to a large scale, which enables a precise and high-resolution spatiotemporal pandemic simulation for entire countries. MEmilio is still under developement but it is available as Open Source and we encourage everyone to make use of it. If you use it, please cite:

- Kühn, Martin Joachim und Abele, Daniel und Kerkmann, David und Korf, Sascha Alexander und Zunker, Henrik und Wendler, Anna Clara und Bicker, Julia und Nguyen, Dang Khoa und Klitz, Margrit und Koslow, Wadim und Siggel, Martin und Kleinert, Jan und Rack, Kathrin und Binder, Sebastian und Plötzke, Lena und Schmieding, René und Lenz, Patrick und Betz, Maximilian Franz und Lutz, Annette und Gerstein, Carlotta und Schmidt, Agatha und Meyer-Hermann, Michael und Basermann, Achim  (2022) MEmilio - a high performance Modular EpideMIcs simuLatIOn software (2022). https://github.com/DLR-SC/memilio, https://elib.dlr.de/192140/.

- Koslow W, Kühn MJ, Binder S, Klitz M, Abele D, et al. (2022) Appropriate relaxation of non-pharmaceutical interventions minimizes the risk of a resurgence in SARS-CoV-2 infections in spite of the Delta variant. PLOS Computational Biology 18(5): e1010054. https://doi.org/10.1371/journal.pcbi.1010054

**Getting started**

MEmilio builds upon different model types, equation-based and agent-based. Furthermore, there are hybrid, graph-ODE-based models. Among the equation-based models, we provide ordinary differential equation (ODE) and integro-differential equation (IDE) based models. In order to provide highly efficient model implementations, MEmilio builds upon a C++ backend for its model and simulation-related content. Data acquisition, plotting, and machine-learnt models are provided in Python.

Details of the C++ implementation of the epidemiological models can be found in the cpp directory (see the [README](cpp/README.md) there). 

Some regularly used data for simulations of a pathogen's spread in Germany, like contact and inter-county mobility, can be found in [data](data/README.md).

In pycode, different MEmilio python packages are defined. Via our [memilio-simulation](pycode/memilio-simulation) package, you can run our C++ backend from python; this package uses [pybind11](https://github.com/pybind/pybind11) to bind our C++ model code. The [memilio-epidata](pycode/memilio-epidata) package provides tools to download and structure important data such as infection or mobility data. More about the python packages can be found in [Python README](pycode/README.rst).

**Documentation**

Each important part of the project described above is described in detail in the README in the corresponding directory. The README contains e.g. configuration and usage instructions for users and developers.

Also, the code is documented with doxygen and instructions on how to obtain it can be found in the docs folder.
The documentation of the code of the main branch can be found at the following URL:

https://dlr-sc.github.io/memilio/documentation/index.html

**Installation, Usage and Requirements**

Each part has different requirements and usage. Detailed instruction can be found in the corresponding READMEs.

**Development**

* [Git workflow and change process](https://github.com/DLR-SC/memilio/wiki/git-workflow)
* [Coding Guidelines](https://github.com/DLR-SC/memilio/wiki/coding-guidelines)
