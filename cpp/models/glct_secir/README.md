# GLCT SECIR model

This model is based on the Generalized Linear Chain Trick (GLCT). 

The GLCT provides the option to use phase-type distributed stay times in the compartments through the use of subcompartments. The Generalized Linear Chain Trick is an extension of the Linear Chain Trick (as the name already suggests). Phase-type distributions are dense in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type distribution of arbitrary precision can be identified.
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The GLCT model can still be described by an ordinary differential equation system.

For the concept see 
- P. J. Hurtado and A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)
- P. J. Hurtado and C. Richards, "Building mean field ODE models using the generalized linear chain trick & Markov chain theory", 2021. (https://doi.org/10.1080/17513758.2021.1912418)

Here, the eight compartments 
- `Susceptible` ($S$), may become exposed at any time
- `Exposed` ($E$), becomes infected after some time
- `InfectedNoSymptoms` ($I_{NS}$), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` ($I_{Sy}$), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` ($I_{Sev}$), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` ($I_{Cr}$), becomes Recovered or Dead after some time
- `Recovered` ($R$)
- `Dead` ($D$)

are used to simulate the spread of the disease. 
It is possible to include phase-type distributed stay times for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.

## Model equations
Below is an overview of the model variables and the model equations are stated. For a simpler description let $\mathcal{Z}=\{{E,I_{NS},I_{Sy},I_{Sev},I_{Cr}}\}$ be the set of the compartments that can be divided in subcompartments.

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Average number of contacts of a person per day. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the susceptible compartments. |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of nonsymptomatically infected people who are not isolated. |
| $\xi_{I_{Sy}}$               | `RiskOfInfectionFromSymptomatic`                | Proportion of infected people with symptoms who are not isolated. |
| $n_{z}$                         |  `Num(...)`  | Number of subcompartments of compartment $z\in\mathcal{Z}$. (...) refers to the C++-name of $z$ as stated above.|
| $\boldsymbol{\alpha_{z}}$                    |  `StartingProbabilities(...)`               | Vector of size $n_{z}$ with the initial probability of starting in any of the subcompartments of compartment $z\in\mathcal{Z}$. The entries should sum to $1$. |
| $\mathbf{A_{z}^{*}}$                    |  `TransitionMatrix(...z)To(...*)`               | Matrix describing the transitions in between of the subcompartments of $z\in\mathcal{Z}$ that describes the transition to the compartment $*$. |

![equations](https://github.com/SciCompMod/memilio/assets/70579874/e1da5e1d-e719-4c16-9f14-45374be7c353)

Note that the notation $\mathbf{z}(t)$ for $z\in\mathcal{Z}$ stands for a vector. If several transitions are possible from a compartment, the vector is split in order to be able to select the stay times until the transitions in each case phase-type distributed. 
For example, the order $\mathbf{I_{\text{NS}}}(t)=[\mathbf{I_{\text{NS}}^{\text{Sy}}}(t),\mathbf{I_{\text{NS}}^{\text{R}}}(t)]^{T}$ is used. Similar holds true for the other compartments $\mathcal{Z}$. 

Implicitly, the matrices $\mathbf{A_{z}^{*}}$ for one $z\in\mathcal{Z}$ are a block of a big matrix $\mathbf{A_{z}}$ corresponding to the whole vector $\mathbf{z}(t)$. As we have no transitions in between the strains defined for different transition probabilities, we would have many zeros in the matrix. The matrix can be defined as

```math
\mathbf{A_{z}}=
\begin{bmatrix}
\mathbf{A_{z}^{*_1}} &  \mathbf{0} \\
\mathbf{0} &  \mathbf{A_{z}^{*_2}}
\end{bmatrix},
```

where ${\*}\_{1}$ is the compartment of the first transition, e.g. $I_{\text{Sy}}$ for $z=I_{\text{NS}}$ and $*_2$ the compartment of the second possible transition, e.g. $R$.
Therefore, we just store the non-zero blocks of the matrix.
Using these parameters, the phase-type distribution that defines the stay time in compartment $z\in\mathcal{Z}$ has the probability density function

$f(x)=\boldsymbol{\alpha_z}^T\hspace{3mu}e^{x\hspace{0.1em}\mathbf{A_z}}\hspace{3mu}(-\mathbf{A_z}\hspace{3mu}\boldsymbol{\Bbb{1}})$ for $x\in\mathbb{R}^{+}$

and the cumulative distribution function

$F(x)=1-\boldsymbol{\alpha_z}^T\hspace{3mu}e^{x\hspace{0.1em}\mathbf{A_z}}\hspace{3mu}\boldsymbol{\Bbb{1}},$

where $e^{x\hspace{0.1em}\mathbf{A_z}}=\sum_{j=0}^{\infty}\frac{(x\hspace{2.5mu}\mathbf{A_z})^j}{j!}$ is the matrix exponential and $\boldsymbol{\Bbb{1}}$ is the vector containing ones of the matching size. Therefore, via changing the vector $\boldsymbol{\alpha_z}$ and the matrices $\mathbf{A_{z}^{*}}$, one can choose the stay time distribution appropriately.

It is important that the sizes of the vectors and matrices match each other and satisfy some other conditions that are checked before a simulation.

The compartment structure with subcompartments is the same as in the LCT-SECIR model. An overview of the model architecture can be found in the [README of the LCT model](../lct_secir/README.md). 
For the GLCT model, some additional transitions are possible and we have more arrows in the model architecture. Below is an example for the Exposed compartment. Note that some Indices are omitted (e.g. $n$ instead of $n_E$) to keep the picture simple.

![tikzGLCTSECIR](https://github.com/user-attachments/assets/fc075b7a-6cd2-4e70-bdd0-a2f4b9f2cf53)
## Examples

A simple example can be found at [GLCT minimal example](../../examples/glct_secir.cpp).
