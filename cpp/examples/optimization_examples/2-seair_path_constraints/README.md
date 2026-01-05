# SEAIR Optimal Control with Path Constraints

This directory contains optimization examples for a **SEAIR epidemiological model** with **time-dependent control interventions** and **path constraints** enforced over the simulation horizon. The examples demonstrate different ways of formulating path constraints in IPOPT and highlight their impact on problem size and solver performance.

## Problem Description

We consider a SEAIR compartmental model with the following **control variables**:

- **SocialDistancing**
- **Quarantined**
- **TestingRate**

These controls are piecewise constant over a fixed number of control intervals and are optimized to balance intervention costs against epidemic outcomes.

The system dynamics are simulated using the `memilio` SEAIR ODE model, and the resulting trajectories are embedded directly into a nonlinear program (NLP) solved with IPOPT.

## Path Constraint of Interest

A key modeling requirement is to **limit the number of infected individuals** during the entire simulation horizon:

Infected(t) ≤ 100,000 for all t ∈ [t₀, t_max]

## Constraint Formulations

### 1. Individual Path Constraints (`seair_pc_individual`)

In this formulation, the path constraint is enforced **at every time step**:

Infected(tₖ) ≤ 100,000 for k = 1, …, N

#### Characteristics
- One constraint per time step
- Total number of constraints scales with the time discretization
- Large constraint vector `m`
- Higher memory usage
- Slower IPOPT convergence due to:
  - Larger KKT system
  - More Lagrange multipliers
  - Increased nonlinear constraint evaluations

This formulation is straightforward and intuitive but does **not scale well** for long horizons or fine time grids.

---

### 2. Global Maximum Path Constraint (`seair_pc_global_max`)

In this formulation, the maximum number of infected individuals over the entire horizon is computed first:

max_t Infected(t) ≤ 100,000

Only **a single scalar constraint** is passed to IPOPT.

#### Characteristics
- Single path constraint (`m` is much smaller)
- Significantly reduced NLP size
- Faster constraint evaluation
- Fewer IPOPT iterations required for convergence
- Better numerical behavior

#### Why this is Better

- IPOPT’s performance is strongly affected by the number of constraints.
- Reducing `m` leads to:
  - Smaller Jacobians
  - Smaller KKT systems
  - Faster linear solves
- The global maximum constraint captures the *same physical requirement* while being computationally more efficient.

➡ **This is the recommended approach for large-scale optimal control problems with path constraints.**

## Objective Function

The objective penalizes intervention intensity over time:

- Social distancing and quarantining are penalized negatively
- Testing is weakly incentivized

The total cost is accumulated over all control intervals.

## Numerical Setup

- Controls are piecewise constant
- Each control interval is subdivided into multiple simulation steps
- SEAIR dynamics are integrated using `memilio::simulate`
- Constraints and objective are evaluated using template-based code to support:
  - double precision
  - automatic differentiation types

## Output

After optimization, the following files are written:

- **State trajectories** over time (CSV)
- **Optimal control profiles** per control interval (CSV)

These can be used for post-processing and visualization.

## Key Takeaway

This example demonstrates that **how path constraints are formulated matters**:

- Enforcing constraints at every time step is simple but inefficient
- Aggregating path constraints via a global maximum:
  - Preserves correctness
  - Improves performance
  - Scales significantly better

This pattern is broadly applicable to optimal control problems beyond epidemiological models.
