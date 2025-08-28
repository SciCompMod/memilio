# ODE MSEIRS4 Model

Implementation of the MSEIRS4 model (maternal immunity + susceptible stages S1-S4 + exposed, infectious, recovered stages with four parallel infection/recovery classes).
This model is designed for Respiratory Syncytial Virus (RSV) and is based on:

- Weber A, Weber M, Milligan P. (2001). *Modeling epidemics caused by respiratory syncytial virus (RSV).* Mathematical Biosciences 172(2): 95–113. `DOI:10.1016/S0025-5564(01)00066-9 <https://doi.org/10.1016/S0025-5564(01)00066-9>`_


## Meaning of indices 1–4 (S1–S4, E1–E4, I1–I4, R1–R4)

- S1: fully susceptible after loss of maternal immunity (highest susceptibility).
- S2: susceptible after first infection (R1 → S2 via waning; reduced susceptibility).
- S3: susceptible after second infection (R2 → S3).
- S4: susceptible after ≥3 infections (R3 → S4 and R4 → S4; lowest susceptibility).

Correspondingly, E_k/I_k/R_k represent exposed/infectious/recovered states of the k-th infection episode. All infectious classes (I1..I4) contribute to transmission equally in the basic formulation of the paper.
