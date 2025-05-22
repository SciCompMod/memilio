GLCT SECIR model
================

This model is based on the Generalized Linear Chain Trick (GLCT).

The GLCT provides the option to use phase-type distributed stay times in the compartments through the use of subcompartments. The Generalized Linear Chain Trick is an extension of the Linear Chain Trick (as the name already suggests). Phase-type distributions are dense in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type distribution of arbitrary precision can be identified. The normal ODE models have (possibly unrealistic) exponentially distributed stay times. The GLCT model can still be described by an ordinary differential equation system.

For the concept see:

- `P. J. Hurtado and A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019 <https://doi.org/10.1007/s00285-019-01412-w>`_
- `P. J. Hurtado and C. Richards, "Building mean field ODE models using the generalized linear chain trick & Markov chain theory", 2021 <https://doi.org/10.1080/17513758.2021.1912418>`_

Here, the eight compartments

- **Susceptible** (:math:`S`), may become exposed at any time
- **Exposed** (:math:`E`), becomes infected after some time
- **InfectedNoSymptoms** (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
- **InfectedSymptoms** (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
- **InfectedSevere** (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
- **InfectedCritical** (:math:`I_{Cr}`), becomes Recovered or Dead after some time
- **Recovered** (:math:`R`)
- **Dead** (:math:`D`)

are used to simulate the spread of the disease. It is possible to include phase-type distributed stay times for the five compartments *Exposed*, *InfectedNoSymptoms*, *InfectedSymptoms*, *InfectedSevere* and *InfectedCritical*.

Model equations
---------------

Below is an overview of the model variables and the model equations are stated. For a simpler description let :math:`\mathcal{Z}=\{E,I_{NS},I_{Sy},I_{Sev},I_{Cr}\}` be the set of the compartments that can be divided into subcompartments.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Average number of contacts of a person per day.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in the Susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of nonsymptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``RiskOfInfectionFromSymptomatic``
     - Proportion of infected people with symptoms who are not isolated.
   * - :math:`n_{z}`
     - ``Num(...)``
     - Number of subcompartments of compartment :math:`z \in \mathcal{Z}`. (...) refers to the C++-name of :math:`z` as stated above.
   * - :math:`\boldsymbol{\alpha_{z}}`
     - ``StartingProbabilities(...)``
     - Vector of size :math:`n_{z}` with the initial probability of starting in any of the subcompartments of compartment :math:`z \in \mathcal{Z}`. The entries should sum up to 1.
   * - :math:`\mathbf{A_{z}^{*}}`
     - ``TransitionMatrix(...z)To(...*)``
     - Matrix describing the transitions in between the subcompartments of :math:`z \in \mathcal{Z}` that describes the transition to the compartment *.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/e1da5e1d-e719-4c16-9f14-45374be7c353
   :alt: equations

Note that the bold notation :math:`\mathbf{z}(t)` for :math:`z \in \mathcal{Z}` stands for a vector. If several transitions are possible from a compartment, the vector is split in order to be able to select the stay times until the transitions in each case phase-type distributed. For example, the order

.. math::

   \mathbf{I_{\text{NS}}}(t) = \begin{bmatrix}
   \mathbf{I_{\text{NS}}^{\text{Sy}}}(t) \\
   \mathbf{I_{\text{NS}}^{\text{R}}}(t)
   \end{bmatrix}

is used. Similar holds true for the other compartments :math:`z \in \mathcal{Z}`.

Implicitly, the matrices :math:`\mathbf{A_{z}^{*}}` for one :math:`z \in \mathcal{Z}` are a block of a matrix :math:`\mathbf{A_{z}}` corresponding to the whole vector :math:`\mathbf{z}(t)`. As we have no transitions in between the strains defined for different transition probabilities, we would have many zeros in the matrix. The matrix can be defined as

.. math::

   \mathbf{A_{z}}=
   \begin{bmatrix}
   \mathbf{A_{z}^{*_1}} &  \mathbf{0} \\
   \mathbf{0} &  \mathbf{A_{z}^{*_2}}
   \end{bmatrix},

where :math:`{*}_{1}` is the compartment of the first transition, e.g., :math:`I_{\text{Sy}}` for :math:`z=I_{\text{NS}}` and :math:`*_{2}` the compartment of the second possible transition, e.g., :math:`R`. Therefore, we just store the non-zero blocks of the matrix. Using these parameters, the phase-type distribution that defines the stay time in compartment :math:`z \in \mathcal{Z}` has the probability density function

.. math::

   f(x)=\boldsymbol{\alpha_z}^T\, e^{x\,\mathbf{A_z}}\, \Bigl(-\mathbf{A_z}\,\boldsymbol{\Bbb{1}}\Bigr)
   \quad \text{for } x\in\mathbb{R}^{+}

and the cumulative distribution function

.. math::

   F(x)=1-\boldsymbol{\alpha_z}^T\, e^{x\,\mathbf{A_z}}\, \boldsymbol{\Bbb{1}},

where

.. math::

   e^{x\,\mathbf{A_z}}=\sum_{j=0}^{\infty}\frac{\bigl(x\,\mathbf{A_z}\bigr)^j}{j!}

is the matrix exponential and :math:`\boldsymbol{\Bbb{1}}` is the vector containing ones of the matching size. Therefore, by changing the vector :math:`\boldsymbol{\alpha_z}` and the matrices :math:`\mathbf{A_{z}^{*}}`, one can choose the stay time distribution appropriately.

It is important that the sizes of the vectors and matrices match each other and satisfy some other conditions that are checked before a simulation.

The compartment structure with subcompartments is the same as in the LCT-SECIR model. An overview of the model architecture can be found in the
`LCT model <lsecir>`_. For the GLCT model, some additional transitions are possible and we have more arrows in the model architecture. Below is an example for the Exposed compartment. Note that some indices are omitted (e.g., :math:`n` instead of :math:`n_E`) to keep the picture simple.

.. image:: https://github.com/user-attachments/assets/fc075b7a-6cd2-4e70-bdd0-a2f4b9f2cf53
   :alt: tikzGLCTSECIR

Examples
--------

A simple example can be found at the
`GLCT minimal example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/glct_secir.cpp>`_.



Overview of the ``glsecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::glsecir