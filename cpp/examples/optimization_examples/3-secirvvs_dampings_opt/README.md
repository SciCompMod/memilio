# SECIRVVS Optimal Control with Damping-Based Interventions

This directory contains optimal control examples for a SECIRVVS epidemiological model using damping factors as control variables. The controls represent non-pharmaceutical interventions (NPIs) that reduce contact intensities in different social settings and are optimized subject to path and terminal constraints.

# Note

This script is intended for testing and prototyping purposes.  
The fully developed version used in **ESID** can be found in:

- `cpp/tools/secirvvs_ESID_dampings`
- `cpp/scripts/optimal-control-applications/4-tools-secirvvs_ESID_dampings`

# Control Variables (Dampings)

Instead of directly controlling rates or compartments, 
we optimize damping factors that scale contact matrices and transmission pathways:

- SchoolClosure
- HomeOffice
- PhysicalDistancingSchool
- PhysicalDistancingWork
- PhysicalDistancingOther

# Interpretation

Each control takes values in [0, 1]
A value of:
- 0 -> no damping (normal contacts)
- 1 -> full damping (maximum reduction)

Controls are piecewise constant over predefined control intervals (e.g., weekly)

These dampings act as continuous policy intensities, 
making the optimization smoother and more realistic than binary on/off interventions.
