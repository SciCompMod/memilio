#!/usr/bin/env python3
"""Benchmark ODE-SECIR mobility through the Python bindings."""
import argparse
import time

try:
    import resource
except ImportError:  # pragma: no cover - only relevant outside POSIX clusters
    resource = None

import numpy as np

import memilio.simulation as mio
import memilio.simulation.osecir as osecir


def make_model(num_groups, seed_scale):
    model = osecir.Model(num_groups)
    p, fact = model.parameters, 1.0 / num_groups
    p.StartDay, p.Seasonality.value = 0.0, 0.0
    p.TestAndTraceCapacity.value = 1e12
    p.ContactPatterns.cont_freq_mat[0].baseline = np.full(
        (num_groups, num_groups), 10.0 * fact)
    for group in range(num_groups):
        age = mio.AgeGroup(group)
        p.TimeExposed[age], p.TimeInfectedNoSymptoms[age] = 3.2, 2.0
        p.TimeInfectedSymptoms[age], p.TimeInfectedSevere[age] = 5.8, 9.5
        p.TimeInfectedCritical[age] = 7.1
        p.TransmissionProbabilityOnContact[age] = 0.1
        p.RelativeTransmissionNoSymptoms[age] = 0.7
        p.RecoveredPerInfectedNoSymptoms[age] = 0.09
        p.RiskOfInfectionFromSymptomatic[age] = 0.25
        p.MaxRiskOfInfectionFromSymptomatic[age] = 0.45
        p.SeverePerInfectedSymptoms[age], p.CriticalPerSevere[age] = 0.2, 0.25
        p.DeathsPerSevere[age], p.DeathsPerCritical[age] = 0.0, 0.3
        for state, value in ((osecir.InfectionState.Exposed, 100.0),
                             (osecir.InfectionState.InfectedNoSymptoms, 50.0),
                             (osecir.InfectionState.InfectedSymptoms, 50.0),
                             (osecir.InfectionState.InfectedSevere, 20.0),
                             (osecir.InfectionState.InfectedCritical, 10.0),
                             (osecir.InfectionState.Recovered, 10.0)):
            model.populations[age, state] = seed_scale * value * fact
        model.populations.set_difference_from_group_total_AgeGroup(
            (age, osecir.InfectionState.Susceptible), 10000.0 * fact)
    model.apply_constraints()
    return model


def mobility(num_groups, degree):
    coefficients = np.full(10 * num_groups, 0.01 / degree if degree else 0.0)
    coefficients.reshape(num_groups, 10)[
        :, int(osecir.InfectionState.Dead)] = 0.0
    return mio.MobilityParameters(coefficients)


def count_steps(tmax, dt):
    steps, current = 0, 0.0
    while current < tmax:
        step = tmax - current if current + dt > tmax else dt
        current += step
        steps += 1
    return steps


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("nodes", type=int)
    parser.add_argument("age_groups", type=int)
    parser.add_argument("tmax", type=float)
    parser.add_argument("dt", type=float)
    parser.add_argument("state_file", nargs="?")
    args = parser.parse_args()
    if (args.nodes < 5 or args.age_groups < 1 or args.tmax <= 0 or args.dt <= 0):
        parser.error("require N >= 5, AGE_GROUPS > 0, TMAX > 0 and DT > 0")
    mio.set_log_level(mio.LogLevel.Critical)

    edges = (2 * args.nodes * (args.nodes - 1) + 2) // 5
    base_degree, extra = divmod(edges, args.nodes)
    all_begin = time.perf_counter()
    model_a = make_model(args.age_groups, 1.0)
    model_b = make_model(args.age_groups, 0.5)
    model_end = time.perf_counter()

    graph = osecir.MobilityGraph()
    for node in range(args.nodes):
        graph.add_node(node, model_a if node % 2 == 0 else model_b, 0.0,
                       args.dt)
    nodes_end = time.perf_counter()
    for node in range(args.nodes):
        graph.get_node(node).property.integrator = mio.RK4IntegratorCore()
    integrator_end = time.perf_counter()
    mobility_by_degree = (mobility(args.age_groups, base_degree),
                          mobility(args.age_groups, base_degree + 1))
    # Match the sorted cyclic-successor topology used by the C++ runner.
    for source in range(args.nodes):
        degree = base_degree + (source < extra)
        params = mobility_by_degree[source < extra]
        wrapped = max(0, source + degree - args.nodes + 1)
        for target in range(wrapped):
            graph.add_edge(source, target, params)
        for target in range(source + 1, min(args.nodes, source + degree + 1)):
            graph.add_edge(source, target, params)
    graph_end = time.perf_counter()

    simulation = osecir.MobilitySimulation(graph, 0.0, args.dt)
    sim_end = time.perf_counter()
    simulation.advance(args.tmax)
    advance_end = time.perf_counter()
    state = np.empty((args.nodes, 10 * args.age_groups))
    for node in range(args.nodes):
        result = simulation.graph.get_node(node).property.result
        state[node] = result.get_last_value()
    extract_end = time.perf_counter()
    state_sum = state.sum(dtype=np.float64)

    if args.state_file:
        state.tofile(args.state_file)
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss if resource else 0
    values = ("python", args.nodes, args.age_groups,
              edges / (args.nodes * (args.nodes - 1)), edges,
              args.tmax, args.dt, count_steps(args.tmax, args.dt), "rk4",
              model_end - all_begin,
              nodes_end - model_end, integrator_end - nodes_end,
              graph_end - integrator_end,
              sim_end - graph_end, sim_end - all_begin,
              advance_end - sim_end, extract_end - advance_end,
              extract_end - all_begin, rss, state_sum)
    print(",".join(str(value) for value in values))


if __name__ == "__main__":
    main()
