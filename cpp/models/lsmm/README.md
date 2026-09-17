# Local stochastic model with endpoint reconstruction

This model implements Pseudocode S3: simulate the aggregate CTMC once, update the
integer history matrix `H` at every event, and sample group endpoints only at the
end of the interval. The transition graph may have branches and cycles.

See [`cpp/examples/lsmm.cpp`](../../examples/lsmm.cpp) for a complete SIR example.
Build and run `lsmm_example` using the usual MEmilio CMake build.

## Exact simulation

- `Model(C, transitions)` describes single-person transitions. Each transition
  has a source index, a target index, and a nonnegative **per-capita** rate
  `double(const State& Z)`. All groups share these rates.
- An optional fourth field lists the state indices used by that rate callback.
  Use `std::vector<Eigen::Index>{1}` for a rate depending only on state 1 and
  `std::vector<Eigen::Index>{}` for a constant rate. Omitting the field (or using
  `std::nullopt`) conservatively assumes dependence on every state. Bare `{}`
  also means `std::nullopt`, so use the explicit empty vector for a constant.
  The list must include every state on which the callback depends. The source
  count needed to form the total transition rate is added automatically.
- Initial integer counts have shape `C × G`, with one group per column. Counts
  must be nonnegative; the total population is limited to `INT_MAX` to keep the
  sampler's integer products in range. Empty groups and entirely empty populations are valid.
- `Simulation(model, initial, rng, t0)` starts with `Z = initial.rowwise().sum()`
  and `H = diag(Z)`. The supplied MEmilio RNG must outlive the simulation.
- `advance(t1)` uses an internal-time next-reaction algorithm. It stores the
  remaining integrated propensity for every channel and scans a vector to select
  the next event. Precomputed dependency lists select the rates to update after
  each event. Per-capita rates are cached until a declared dependency changes,
  while changed source counts always update the total rates. Calls to `advance`
  preserve both the rate caches and the residual clocks. An event at the end
  time is included. Clock accounting and event selection still scan all channels.
- An integer categorical draw from the source row of `H` determines the initial
  state of the affected individual. That column is updated immediately. No
  transition counts or marked event sequence are used for reconstruction.
- `sample_endpoints()` samples the endpoints conditional on `H`. It uses
  sequential univariate hypergeometric draws to form the multivariate draw for
  each group and initial state. The last group receives the residual counts.
  The same sampler is available as `reconstruct(initial, H, rng)`.

Rates are autonomous functions of the aggregate state; the callback has no time
argument. Do not mutate captured rate parameters during an interval. Prescribed
piecewise-constant changes and discrete transport are handled by ending the
interval, sampling endpoints, applying the change/transport, and constructing a
new simulation. General continuously time-dependent hazards and delays are not
implemented. Neither births nor simultaneous changes of multiple individuals
are allowed.

The SIR example captures the conserved population size once, declares infection
as dependent on `I`, and marks the per-capita recovery rate as constant. Both
total rates, `beta*S*I/N` and `gamma*I`, still need updating after infection or
recovery. Caching `gamma` does not make the total recovery rate constant.

`H(current, initial)` has row sums equal to the current aggregate and column sums
equal to the initial aggregate. Group totals are conserved exactly. Repeated
endpoint draws from the same `H` are alternative conditional samples; they do
not constitute a consistent group-resolved time series.

Core simulation storage is `O(C² + CG + R)`, plus the supplied state dependencies
and the precomputed lists of affected channels (at most `R²` channel pairs).
There is no group loop during the aggregate event update.
The endpoint sampler uses mode-centered hypergeometric inversion with adjacent
probability ratios. Relative masses are normalized in long double; geometric
tail bounds stop at machine precision. No factorials, growing event sequence,
or per-individual reconstruction loop is needed. The number of scalar draws is
`O(GC²)`; each draw's work depends on the width of its distribution, so scalar
sampling should not be assumed to have constant cost.

## Optional approximation

`simulate_tau_leaping(model, initial, rng, t0, t1, dt)` is a separate **approximate
bounded binomial tau-leaping** implementation. It freezes the rates at the start
of each step, draws departures from each history cell with probability
`1 - exp(-total_exit_rate * dt)`, and distributes departures over outgoing
channels using their relative rates. All updates use the old history matrix, so
each individual changes state at most once per step. Nonnegative integer counts,
history margins, and group totals are preserved, even with branches and cycles.

The result contains `state`, `history`, and sampled `endpoints`. This is not the
exact CTMC of S3: it omits multiple successive transitions within one step and
freezes state-dependent rates. Check time-step convergence for the intended
model. Endpoint reconstruction is exact conditional on the approximate history;
it does not remove the time-stepping error. No speedup is guaranteed.
