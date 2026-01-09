# Optimization Under Uncertainty

This tool extends the SECIRVVS ESID damping optimization to account for uncertain model parameters.

It is based on the example 4-tools-secirvvs-ESID-dampings example.

# Uncertainty model

The optimization accounts for uncertainty in:

- Transmission probability. 
Sampled from the interval [0.07, 0.12]

- Initial infected population states.
Perturbations applied to the infected compartments

Each optimization step runs the epidemic model multiple times using different parameter samples.

# Optimization objective

The control problem:

- Minimizes the average objective value across all sampled runs
- Enforces constraints for every run, not just on average

A control schedule is accepted only if all realizations satisfy the constraints.

# Interpretation

The resulting control strategy is:

- Robust to parameter uncertainty
- Conservative with respect to constraint violations
- Optimized for expected performance rather than a single scenario
