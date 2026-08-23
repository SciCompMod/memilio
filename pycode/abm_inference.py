"""Shared ABM/BayesFlow inference building blocks
"""
import os

if "KERAS_BACKEND" not in os.environ:
    os.environ["KERAS_BACKEND"] = "jax"

from collections.abc import Collection

import matplotlib.pyplot as plt
import numpy as np
import keras
import bayesflow as bf
from scipy import stats

from memilio.simulation.abm import TestingBudget
from abm_batch import run_batch_designs

# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────
TOTAL_POPULATION = 100_000

SIM_DAYS = 30                    # true simulated horizon (must match forward_pass.cpp's hardcoded
                                 # tmax -- the two aren't linked, so keep them in sync by hand).
                                 # Governs n_events()/make_design()'s budget accounting and how much
                                 # of the record log is in-range at all; NOT how much of it is used
                                 # as the observable -- see OBSERVED_DAYS.
REPORTING_DELAY_DAYS = 1         # must match TestingBudget's default reporting_delay
                                 # (pcr_surveillance.h); make_design() never overrides it. Update
                                 # this (and re-derive OBSERVED_DAYS's margin) if that ever changes.
OBSERVED_DAYS = SIM_DAYS - REPORTING_DELAY_DAYS  # how many of SIM_DAYS' days go into
                                 # histogram_ct/histogram_bin. Kept short of SIM_DAYS so every
                                 # included day's tests have had time to resolve: a test taken on
                                 # day d only appears in forward_pass()'s output once
                                 # resolve_pending() fires at d + reporting_delay, which for the
                                 # last REPORTING_DELAY_DAYS simulated day(s) falls at or past tmax
                                 # and so never happens -- those days would otherwise look
                                 # artificially, uniformly empty in every simulation regardless of
                                 # theta. Verified empirically (day 28 resolves at exactly ratio
                                 # 1.0, day 29 at exactly 0.0 -- a deterministic boundary, not a
                                 # probabilistic edge effect) that this margin needs no slack beyond
                                 # reporting_delay itself.
SURVEY, DIAGNOSTIC = 0, 1        # source tags from forward_pass: 0 = active random surveillance,
                                 # 1 = symptomatic diagnostic care-seeking
N_CT_BINS = 41                   # CT 0..40

# Which test stream(s) to build histogram_ct/histogram_bin from -- pass one of these as
# _assemble()'s/make_val_at_design()'s `sources` to compare survey-only (the default so far),
# diagnostic-only, or both pooled. Diagnostic care-seeking already runs in every simulation
# (it suppresses the epidemic via isolation-on-positive) but its own test results weren't
# previously used as an observable -- these make that comparison a one-argument change.
SOURCES_SURVEY = frozenset({SURVEY})
SOURCES_DIAGNOSTIC = frozenset({DIAGNOSTIC})
SOURCES_POOLED = frozenset({SURVEY, DIAGNOSTIC})

PARAM_KEYS = ["beta", "t_exposed", "time_presymptomatic", "time_asymptomatic_recovery", "symptom_prob",
             "quarantine_compliance"]
PAR_NAMES = [r"$\beta$", r"$t_\mathrm{exp}$", r"$t_\mathrm{presym}$", r"$t_\mathrm{asym}$", r"$p_\mathrm{sym}$",
            r"$p_\mathrm{quar}$"]

# theta priors -- full support (lognormal for positive params, logit-normal for the
# probabilities), exactly matching joint_inference.ipynb's priors.
PRIOR_BETA_S,   PRIOR_BETA_SCALE  = 0.5,  1.0   # median 1.0
PRIOR_TEXP_S,   PRIOR_TEXP_SCALE  = 0.4,  5.0   # median 5.0  (~95% in [2.3, 11])
PRIOR_TPRE_S,   PRIOR_TPRE_SCALE  = 0.5,  1.5   # median 1.5  (~95% in [0.56, 4.0])
PRIOR_TASYM_S,  PRIOR_TASYM_SCALE = 0.25, 8.0   # median 8.0  (~95% in [4.9, 13])
PRIOR_PSYM_MU,  PRIOR_PSYM_SIGMA  = np.log(0.6 / 0.4), 0.5  # logit-normal, median 0.6
# Was a fixed city-level constant (QUARANTINE_COMPLIANCE = 0.6); now inferred like every other
# theta. Same logit-normal shape/spread as symptom_prob, centered at the old fixed value so the
# prior's mass is where the model used to assume it was, but now with real (0,1) uncertainty
# around it -- (~95% in [0.35, 0.80]).
PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA = np.log(0.6 / 0.4), 0.5  # logit-normal, median 0.6

NUM_POSTERIOR_SAMPLES = 500


# ─────────────────────────────────────────────────────────────────────────────
# Data-generating process
# ─────────────────────────────────────────────────────────────────────────────
def _daily_test_totals(out, sources, sim_days):
    """(sim_days,) total number of tests taken per day (positives + negatives), restricted to
    `sources`. Used as records_to_histogram()'s rate denominator."""
    totals = np.zeros(sim_days)
    pos = out["positives"]
    if pos.shape[0]:
        m = np.isin(pos[:, 5], list(sources))
        days = np.round(pos[m, 0]).astype(int)
        valid = (days >= 0) & (days < sim_days)
        np.add.at(totals, days[valid], 1)
    neg = out["negatives"]
    if neg.shape[0]:
        neg_days = np.round(neg[:, 0]).astype(int)
        m = np.isin(neg[:, 1], list(sources))
        valid = m & (neg_days >= 0) & (neg_days < sim_days)
        np.add.at(totals, neg_days[valid], neg[valid, 4])
    return totals


def records_to_histogram(out, sources: Collection[int] = frozenset({SURVEY}), sim_days=OBSERVED_DAYS,
                         n_bins=N_CT_BINS, scale_by_tests=False):
    """forward_pass() dict -> (sim_days, n_bins) city-wide daily CT histogram of positives.
    `sources` is a set of stream tags to pool: {0}=survey, {1}=diagnostic, {0, 1}=both.
    Positives are out['positives'] = [day, person, age, loc, ct, source].

    `sim_days` defaults to OBSERVED_DAYS, not the true simulated horizon SIM_DAYS: the model
    actually runs for SIM_DAYS days regardless (this only slices what's read out of the record
    log), but the last few simulated days' tests haven't had time to resolve within the
    simulation's own horizon and would otherwise show up as spuriously, uniformly empty -- see
    OBSERVED_DAYS. Records from day >= sim_days are dropped outright, not merely clipped.

    If `scale_by_tests`, each day's row is divided by that day's total number of tests taken
    (positives + negatives, restricted to the same `sources`) -- turning raw positive COUNTS into
    a per-bin test-POSITIVITY RATE. This implicitly folds negative results back into the
    observable without a separate negative-count channel: a day with many negatives dilutes the
    rate even though the positive count alone is unchanged.
    histogram_bin (city.sum(axis=-1) in _assemble()) then becomes the classic daily
    test-positivity-rate curve. Days with zero tests get an all-zero row (as they already did
    for raw counts, so this isn't a new degenerate case).

    The denominator itself isn't kept as a separate feature: log_budget/frequency already pin
    down the nominal daily survey test count exactly, and empirically the realized count tracks
    it almost everywhere in this model's regime (isolation essentially never depletes the
    eligible survey frame below the nominal sample size), so a separate count channel would
    mostly duplicate what's already implicit in the existing design conditions."""
    city = np.zeros((sim_days, n_bins))
    pos = out["positives"]
    if pos.shape[0]:
        m = np.isin(pos[:, 5], list(sources))
        days = np.round(pos[m, 0]).astype(int)
        cts = np.clip(np.round(pos[m, 4]).astype(int), 0, n_bins - 1)
        valid = (days >= 0) & (days < sim_days)
        np.add.at(city, (days[valid], cts[valid]), 1)

    if not scale_by_tests:
        return city

    totals = _daily_test_totals(out, sources, sim_days)
    return np.divide(city, totals[:, None], out=np.zeros_like(city), where=totals[:, None] > 0)


def prior_theta(rng):
    """One joint draw over all PARAM_KEYS; the returned dict IS a forward_pass params dict."""
    z_sym = rng.normal(PRIOR_PSYM_MU, PRIOR_PSYM_SIGMA)
    z_qcomp = rng.normal(PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA)
    return {
        "beta":                       rng.lognormal(np.log(PRIOR_BETA_SCALE),  PRIOR_BETA_S),
        "t_exposed":                  rng.lognormal(np.log(PRIOR_TEXP_SCALE),  PRIOR_TEXP_S),
        "time_presymptomatic":        rng.lognormal(np.log(PRIOR_TPRE_SCALE),  PRIOR_TPRE_S),
        "time_asymptomatic_recovery": rng.lognormal(np.log(PRIOR_TASYM_SCALE), PRIOR_TASYM_S),
        "symptom_prob":               1.0 / (1.0 + np.exp(-z_sym)),  # expit(N(PSYM_MU, PSYM_SIGMA))
        "quarantine_compliance":      1.0 / (1.0 + np.exp(-z_qcomp)),  # expit(N(QCOMP_MU, QCOMP_SIGMA))
    }


def _logitnorm_logpdf(p, mu, sigma):
    """log-density of a logit-normal at p in (0,1); -inf outside."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, -np.inf)
    m = (p > 0.0) & (p < 1.0)
    z = np.log(p[m] / (1.0 - p[m]))
    out[m] = stats.norm.logpdf(z, loc=mu, scale=sigma) - np.log(p[m]) - np.log1p(-p[m])
    return out


def prior_log_prob(**params):
    """log p(theta) under prior_theta()'s exact independent sampling distributions. Each factor
    is -inf outside its support, which estimate_mutual_information() masks."""
    return (stats.lognorm.logpdf(params["beta"], s=PRIOR_BETA_S,  scale=PRIOR_BETA_SCALE)
            + stats.lognorm.logpdf(params["t_exposed"], s=PRIOR_TEXP_S,  scale=PRIOR_TEXP_SCALE)
            + stats.lognorm.logpdf(params["time_presymptomatic"], s=PRIOR_TPRE_S,  scale=PRIOR_TPRE_SCALE)
            + stats.lognorm.logpdf(params["time_asymptomatic_recovery"], s=PRIOR_TASYM_S, scale=PRIOR_TASYM_SCALE)
            + _logitnorm_logpdf(params["symptom_prob"], PRIOR_PSYM_MU, PRIOR_PSYM_SIGMA)
            + _logitnorm_logpdf(params["quarantine_compliance"], PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA))


def prior_mean(n_mc=200_000, rng=None):
    """Analytic/Monte-Carlo mean of prior_theta() per param. The 4 lognormal params have a
    closed form (scale * exp(s^2/2)); the two logit-normal params' means have none, so they're
    Monte-Carlo estimated."""
    means = {
        "beta":                       PRIOR_BETA_SCALE  * np.exp(PRIOR_BETA_S**2 / 2),
        "t_exposed":                  PRIOR_TEXP_SCALE  * np.exp(PRIOR_TEXP_S**2 / 2),
        "time_presymptomatic":        PRIOR_TPRE_SCALE  * np.exp(PRIOR_TPRE_S**2 / 2),
        "time_asymptomatic_recovery": PRIOR_TASYM_SCALE * np.exp(PRIOR_TASYM_S**2 / 2),
    }
    rng = rng or np.random.default_rng()
    z = rng.normal(PRIOR_PSYM_MU, PRIOR_PSYM_SIGMA, n_mc)
    means["symptom_prob"] = float(np.mean(1.0 / (1.0 + np.exp(-z))))
    z = rng.normal(PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA, n_mc)
    means["quarantine_compliance"] = float(np.mean(1.0 / (1.0 + np.exp(-z))))
    return means


def n_events(period, sim_days=SIM_DAYS):
    """Number of surveillance events (testing days 0, period, 2*period, ... < sim_days)."""
    return (sim_days - 1) // period + 1


def make_design(total_budget, period):
    """(total_budget, period) -> a TestingBudget. `total_budget` is the cumulative fraction of the
    population tested (in expectation) over SIM_DAYS days -- a Python/experimental-design-level
    concept the C++ model itself has no notion of: TestingBudget's own event_budget_fraction is a
    plain per-event fraction, independent of cadence. This is where the two meet: total_budget is
    spread evenly over however many events actually occur (n_events(period)) to get there."""
    event_budget_fraction = total_budget / n_events(period)
    return TestingBudget(event_budget_fraction=float(event_budget_fraction), test_period_days=int(period))


def total_budget_for_rate(rate_per_100k_per_day, sim_days=SIM_DAYS):
    """A real-world testing-rate benchmark (tests / 100,000 people / day -- the unit surveillance
    literature is normally reported in) -> this module's `total_budget` unit (cumulative survey
    tests per person over sim_days, what make_design()/DESIGN_TOTAL_LO/HI/TOTAL_BUDGET_GRID
    actually take). Independent of period/cadence -- total_budget is already the cumulative amount
    over the whole window, this just re-expresses a daily per-capita rate as that cumulative total.
    Do unit conversions here, once, rather than mentally re-deriving them at every call site."""
    return rate_per_100k_per_day / 100_000 * sim_days


def rate_per_100k_per_day(total_budget, sim_days=SIM_DAYS):
    """Inverse of total_budget_for_rate(): this module's total_budget -> tests / 100k people / day."""
    return total_budget / sim_days * 100_000


def _assemble(draws, totals, periods, outs, sources=SOURCES_SURVEY, scale_by_tests=False):
    """`sources` selects which test stream(s) histogram_ct/histogram_bin are built from -- see
    SOURCES_SURVEY/SOURCES_DIAGNOSTIC/SOURCES_POOLED. Both diagnostic and survey testing always
    run inside the simulation regardless (diagnostic care-seeking suppresses the epidemic via
    isolation-on-positive either way); this only controls which stream's *results* become the
    observable the network is trained/evaluated on. `scale_by_tests`: see records_to_histogram().

    No separate test-count channel is added here: log_budget/frequency (below) already pin down
    the *nominal* survey test count on cadence days exactly (event_budget_fraction * population),
    and empirically the *realized* count matches it almost everywhere --
    isolation essentially never depletes the eligible survey frame below the nominal sample size,
    even at compliance=1.0. So an explicit count channel would mostly just make the network
    re-derive what log_budget/frequency + the day index already determine (and for Diagnostic,
    the count is ~95% redundant with the positive count itself, since only symptomatic-infected
    Persons seek diagnostic tests in the first place)."""
    city = np.stack([records_to_histogram(o, sources, scale_by_tests=scale_by_tests) for o in outs])
    data = {k: np.array([d[k] for d in draws], dtype=np.float64).reshape(-1, 1) for k in PARAM_KEYS}
    data.update({
        # design conditions, already in network-friendly form:
        "log_budget": np.log(np.asarray(totals, dtype=np.float64)).reshape(-1, 1),
        "frequency":  (1.0 / np.asarray(periods, dtype=np.float64)).reshape(-1, 1),
        "histogram_ct":  city,
        "histogram_bin": city.sum(axis=-1),
        "n_ever_infected": np.array([o["n_ever_infected"][0, 0] for o in outs], dtype=np.float64),
    })
    return data


def _assemble_dual_stream(draws, totals, periods, outs):
    """Like _assemble(), but keeps Survey and Diagnostic separate -- histogram_ct_{survey,
    diagnostic} and daily_test_total_{survey,diagnostic} -- instead of pre-combining into one
    `sources`/`scale_by_tests` choice. Used by simulation_cache.py so a cached dataset can be
    reprocessed into ANY sources/scale_by_tests combination later without resimulating: summing
    the two streams reproduces SOURCES_POOLED exactly (verified empirically -- survey + diagnostic
    == pooled, bit-for-bit, for both the histogram and the daily test totals, since
    records_to_histogram()/_daily_test_totals() both accumulate additively over disjoint per-
    source subsets of the same records), and dividing by the (optionally stream-summed) totals
    reproduces scale_by_tests=True's rate exactly. All arrays float32 -- this is a storage-facing
    assembly, and the adapter downcasts to float32 for training anyway, so float64 here would only
    double the on-disk footprint for no benefit."""
    city_survey = np.stack([records_to_histogram(o, SOURCES_SURVEY) for o in outs]).astype(np.float32)
    city_diag   = np.stack([records_to_histogram(o, SOURCES_DIAGNOSTIC) for o in outs]).astype(np.float32)
    tot_survey  = np.stack([_daily_test_totals(o, SOURCES_SURVEY, OBSERVED_DAYS) for o in outs]).astype(np.float32)
    tot_diag    = np.stack([_daily_test_totals(o, SOURCES_DIAGNOSTIC, OBSERVED_DAYS) for o in outs]).astype(np.float32)
    data = {k: np.array([d[k] for d in draws], dtype=np.float32) for k in PARAM_KEYS}
    data.update({
        "total_budget": np.asarray(totals, dtype=np.float32),
        "period": np.asarray(periods, dtype=np.float32),
        "n_ever_infected": np.array([o["n_ever_infected"][0, 0] for o in outs], dtype=np.float32),
        "histogram_ct_survey": city_survey,
        "daily_test_total_survey": tot_survey,
        "histogram_ct_diagnostic": city_diag,
        "daily_test_total_diagnostic": tot_diag,
    })
    return data


def make_val_at_design(n, population, total_budget, period, rng, show_progress=False, draws=None,
                       sources=SOURCES_SURVEY, scale_by_tests=False):
    """n prior draws (or the given `draws`, e.g. a repeated fixed theta) at a single fixed
    design -- used as "the validation set" in run_design_sweep.py; the `draws` override also
    supports building a training set for a single-design (non design-amortized) workflow.
    `sources`, `scale_by_tests`: see _assemble()."""
    if draws is None:
        draws = [prior_theta(rng) for _ in range(n)]
    designs = [make_design(total_budget, period)] * n
    outs = run_batch_designs(draws, designs, population, show_progress=show_progress, max_workers=os.cpu_count())
    return _assemble(draws, np.full(n, total_budget), np.full(n, period), outs, sources=sources,
                     scale_by_tests=scale_by_tests)


# ─────────────────────────────────────────────────────────────────────────────
# Model
# ─────────────────────────────────────────────────────────────────────────────
@keras.saving.register_keras_serializable(package="abm_inference")
class GRU(bf.networks.SummaryNetwork):
    def __init__(self, hidden_dim=128, summary_dim=16, dropout=0.2, recurrent_dropout=0.1, **kwargs):
        super().__init__(**kwargs)
        self.gru1 = keras.layers.GRU(hidden_dim, return_sequences=True,
                                     dropout=dropout, recurrent_dropout=recurrent_dropout)
        self.gru2 = keras.layers.GRU(hidden_dim, dropout=dropout, recurrent_dropout=recurrent_dropout)
        self.summary_stats = keras.layers.Dense(summary_dim)

    def call(self, time_series, **kwargs):
        training = kwargs.get("stage") == "training"
        x = self.gru1(time_series, training=training)
        x = self.gru2(x, training=training)
        return self.summary_stats(x)


def build_workflow(summary_key, design_conditioned=True, checkpoint_dir=None, simulator=None):
    """Infers all 5 params (PARAM_KEYS) from a histogram summary. If `design_conditioned`, also
    takes the design (log_budget, frequency) as a direct condition, so one trained network
    handles any (total_budget, period) -- set False for a single-design workflow (no design
    varies, so conditioning on it is pure wasted capacity).

    `simulator`, if given, is a bayesflow Simulator attached to the workflow so fit_online() can
    call it live -- see run_design_sweep_online.py. Unused for fit_offline(); None (the default)
    is fine unless you're building an online-trained workflow.

    If `checkpoint_dir` is given, the best (lowest val_loss) model seen during fit_offline()/
    fit_online() is written there as f'{summary_key}.keras' -- without this, BasicWorkflow never
    persists anything and the trained network only exists in process memory, gone the moment the
    run ends (deleting the results directory doesn't cause this; there was simply never anything
    else to keep). `checkpoint_dir` must already exist -- ModelCheckpoint doesn't create it.

    Full-model (not weights-only) on purpose: the Standardization layer's running mean/var/count
    (bayesflow.networks.helpers.standardization) are ordinary non-trainable Keras weights, but a
    weights-only (.weights.h5) restore into a freshly-built workflow was verified (empirically,
    across separate processes) to silently leave `count` at whatever the new instance's own
    build pass produced instead of the checkpointed value -- e.g. 24 instead of 120 -- which then
    corrupts every downstream standardized input. A full .keras restore replaces the whole
    approximator object wholesale and was verified to reproduce log_prob() bit-for-bit.

    To reload: build_workflow() the same summary_key/design_conditioned again (no throwaway fit
    needed first, unlike weights-only restore) and call the new workflow's .load_approximator()
    (path defaults to checkpoint_dir/checkpoint_name from construction, or pass one explicitly).
    Requires GRU (below) to stay registered via @keras.saving.register_keras_serializable --
    that's what lets a full-model checkpoint reference our custom summary network by name."""
    adapter = (
        bf.adapters.Adapter()
        .convert_dtype("float64", "float32")
        .as_time_series(summary_key)
        .concatenate(PARAM_KEYS, into="inference_variables")
    )
    if design_conditioned:
        adapter = adapter.concatenate(["log_budget", "frequency"], into="inference_conditions")
    adapter = (
        adapter
        .rename(summary_key, "summary_variables")
        .log(["inference_variables", "summary_variables"], p1=True)
    )
    standardize = ["inference_variables", "inference_conditions"] if design_conditioned else ["inference_variables"]
    return bf.BasicWorkflow(
        simulator=simulator,
        adapter=adapter,
        inference_network=bf.networks.CouplingFlow(),
        summary_network=GRU(),
        standardize=standardize,
        checkpoint_filepath=str(checkpoint_dir) if checkpoint_dir is not None else None,
        checkpoint_name=summary_key,
        save_weights_only=False,
        save_best_only=True,
    )


def estimate_mutual_information(workflow, val, x_key, num_posterior_samples=NUM_POSTERIOR_SAMPLES,
                                val_batch=32, cond_keys=("log_budget", "frequency")):
    """I(theta;X|d) for a validation set at a fixed design; the design (log_budget, frequency,
    carried in `val`) is passed as a condition to sample() and log_prob(). Pass cond_keys=()
    for a non-design-conditioned workflow."""
    n_val = len(val["beta"])
    kls = np.empty(n_val)
    n_invalid, n_total = 0, 0
    for start in range(0, n_val, val_batch):
        idx = slice(start, min(start + val_batch, n_val))
        cond = {x_key: val[x_key][idx], **{k: val[k][idx] for k in cond_keys}}
        post = workflow.sample(conditions=cond, num_samples=num_posterior_samples)
        samples = {k: post[k].reshape(-1, 1) for k in PARAM_KEYS}
        b = val["beta"][idx].shape[0]
        tiled = {x_key: np.repeat(val[x_key][idx], num_posterior_samples, axis=0)}
        for k in cond_keys:
            tiled[k] = np.repeat(val[k][idx], num_posterior_samples, axis=0)
        log_q = workflow.log_prob(data={**tiled, **samples})
        log_q = log_q.reshape(b, num_posterior_samples)
        log_p = prior_log_prob(**samples).reshape(b, num_posterior_samples)
        diff = np.where(np.isfinite(log_q) & np.isfinite(log_p), log_q - log_p, np.nan)
        n_invalid += np.isnan(diff).sum()
        n_total += diff.size
        kls[idx] = np.nanmean(diff, axis=1)
    return kls, n_invalid / n_total


def sample_posterior(workflow, val, x_key, num_posterior_samples=NUM_POSTERIOR_SAMPLES,
                     val_batch=32, cond_keys=("log_budget", "frequency")):
    """Posterior draws for every row of `val`, batched like estimate_mutual_information().
    Returns {param: (n_val, num_posterior_samples, 1) for param in PARAM_KEYS} --
    workflow.sample()'s own shape convention, trailing singleton axis included. Keep it: the
    bf.diagnostics.plots functions (recovery/calibration_ecdf/pairs_posterior) use that axis to
    tell "one variable's full sample array" apart from "needs splitting into separate variables"
    (see bayesflow.utils.dict_utils.split_arrays) -- squeezing it away silently explodes into
    num_posterior_samples bogus "variables". Pass cond_keys=() for a non-design-conditioned
    workflow."""
    n_val = len(val["beta"])
    out = {k: np.empty((n_val, num_posterior_samples, 1)) for k in PARAM_KEYS}
    for start in range(0, n_val, val_batch):
        idx = slice(start, min(start + val_batch, n_val))
        cond = {x_key: val[x_key][idx], **{k: val[k][idx] for k in cond_keys}}
        post = workflow.sample(conditions=cond, num_samples=num_posterior_samples)
        for k in PARAM_KEYS:
            out[k][idx] = post[k]
    return out


def local_jacobian(workflow, val, idx, x_key, cond_keys=("log_budget", "frequency"), z=None):
    """Exact local Jacobian d(theta)/dz of the trained flow's generative map (base sample z ->
    raw theta, in PARAM_KEYS order) at one validation instance, linearized at `z` (default: the
    base distribution's mode, z=0).

    Autodiffs (via jax) through the coupling flow's inverse map and the standardizer, both plain
    differentiable keras/jax ops; the adapter's log1p transform runs on concrete numpy arrays (not
    jax-traceable), so its Jacobian (diagonal: d expm1(u)/du = exp(u)) is applied by hand via the
    chain rule to land back in raw theta units.

    Returns (J_raw, theta_at_z): J_raw is (len(PARAM_KEYS), len(PARAM_KEYS)).
    """
    import jax

    n_vars = len(PARAM_KEYS)
    appr = workflow.approximator
    cond, _, _ = appr._prepare_conditions(
        {x_key: val[x_key][idx:idx + 1], **{k: val[k][idx:idx + 1] for k in cond_keys}})

    def theta_adapted(z_vec):
        x = appr.inference_network(z_vec[None, :], conditions=cond, inverse=True, training=False)
        return appr.standardizer.maybe_standardize(
            x, key="inference_variables", stage="inference", forward=False)[0]

    if z is None:
        z = jax.numpy.zeros(n_vars, dtype="float32")  # base-distribution mode
    J_std = np.asarray(jax.jacobian(theta_adapted)(z))
    theta_log1p = np.asarray(theta_adapted(z))
    theta_at_z = np.expm1(theta_log1p)
    J_raw = np.diag(np.exp(theta_log1p)) @ J_std  # chain rule: d expm1(u)/du = exp(u)
    return J_raw, theta_at_z


def local_jacobian_svd(workflow, val, idx, x_key, cond_keys=("log_budget", "frequency"), z=None):
    """SVD of local_jacobian(), in RELATIVE (log-parameter) units, NOT raw physical units --
    d(log theta)/dz = (1/theta) * d(theta)/dz.

    This matters and isn't optional: SVD in raw units isn't scale-invariant (verified: rescaling
    one parameter's units, e.g. days -> hours, changes which direction the SVD calls "sloppiest"
    even though nothing about the actual identifiability changed -- raw singular values mix
    different units (days, a probability, a rate) into one Euclidean norm, so they can't be
    legitimately compared across parameters, or across sweep conditions like different designs or
    compliance levels). Relative units are dimensionless fractional sensitivities (0.2 = "20%
    local uncertainty along this direction"), genuinely comparable either way -- the standard fix
    for exactly this issue in the sloppy-models literature (Gutenkunst et al.).

    Since z ~ Normal(0, I) locally, the local posterior covariance in LOG theta is J_rel @ J_rel.T
    to leading order. So SVD(J_rel): columns of U are directions in log-theta space (PARAM_KEYS
    order) -- e.g. a component pattern near (1, -1, 0, 0, 0) reads as "log(beta) - log(t_exposed)
    is unconstrained", i.e. beta/t_exposed =~ const, a genuine Laurent-monomial confounding
    direction -- and each singular value is the local FRACTIONAL posterior std along that
    direction: large = locally confounded/unidentified, small = locally well-constrained. Compare
    U across several `idx` (or several `z`) to see whether the confounded direction is fixed (a
    single global relation) or rotates (genuinely non-monomial structure).

    Returns (singular_values, U, theta_at_z) -- theta_at_z is still in raw units (the center
    point), only S/U are relative. Caveat: symptom_prob is (0,1)-bounded, not a genuinely
    unbounded positive quantity, so 1/theta as its local scale is an approximation (the exact
    scale-free choice would be its logit) -- matches how the flow itself represents it internally
    (log1p uniformly, no logit special-casing), not a new inconsistency introduced here.
    """
    J_raw, theta_at_z = local_jacobian(workflow, val, idx, x_key, cond_keys, z)
    J_rel = np.diag(1.0 / theta_at_z) @ J_raw
    U, S, _ = np.linalg.svd(J_rel)
    return S, U, theta_at_z


def shared_bounds(val, posts, keys=PARAM_KEYS, pad=0.05, pct=(1, 99)):
    """Per-parameter (lo, hi) axis bounds spanning ground truth (`val`) AND every workflow's
    posterior samples in `posts` (dict of label -> sample_posterior() output) -- so multiple
    workflows' plots share one scale instead of each auto-ranging to its own data.

    Uses percentiles (`pct`), not literal min/max: a normalizing flow's posterior for a weakly-
    identified parameter routinely throws a handful of wild-tail draws, and literal min/max would
    stretch the axis out to fit those few outliers, leaving the actual dense scatter squeezed
    into a speck in the middle -- exactly the "zoomed out, can't see anything" failure. `pad` is
    an extra fractional margin on top of the percentile range."""
    bounds = {}
    for k in keys:
        pooled = np.concatenate([val[k].ravel()] + [post[k].ravel() for post in posts.values()])
        lo, hi = np.percentile(pooled, pct)
        margin = pad * (hi - lo) if hi > lo else 1.0
        bounds[k] = (lo - margin, hi + margin)
    return bounds


def apply_recovery_bounds(fig, bounds, keys=PARAM_KEYS):
    """recovery()'s figure has exactly one axes per variable, in `keys` order; ground truth and
    estimate are the same quantity, so both axes take the same bound."""
    for ax, k in zip(fig.axes, keys):
        ax.set_xlim(bounds[k]); ax.set_ylim(bounds[k])


def apply_pairgrid_bounds(g, bounds, keys=PARAM_KEYS):
    """pairs_posterior()'s PairGrid.axes is a clean (len(keys), len(keys)) grid: column j is
    always variable keys[j] (x-axis), row i is always keys[i] (y-axis, off-diagonal only --
    diagonal panels are univariate density plots, y is a density not the variable)."""
    n = len(keys)
    for col, k in enumerate(keys):
        for row in range(n):
            g.axes[row, col].set_xlim(bounds[k])
    for row, k in enumerate(keys):
        for col in range(n):
            if row != col:
                g.axes[row, col].set_ylim(bounds[k])


def plot_local_identifiability(workflow, post, val, x_key, indices, label,
                               var_pair=("beta", "t_exposed"), bounds=None,
                               cond_keys=("log_budget", "frequency")):
    """Empirical posterior scatter (from sample_posterior()) vs. the flow's exact local 1-sigma
    linearization ellipse, for a handful of validation instances -- shows at a glance whether/how
    the locally sloppy direction rotates, projected onto the 2 parameters in `var_pair`.

    The ellipse is the MARGINAL of the full len(PARAM_KEYS)-dim local covariance J @ J.T onto
    `var_pair`'s 2 coordinates (the corresponding 2x2 submatrix, eigendecomposed) -- valid because
    marginalizing a (locally linearized) Gaussian onto a coordinate subset is exactly that
    submatrix, no need to re-run the flow per pair."""
    i, j = PARAM_KEYS.index(var_pair[0]), PARAM_KEYS.index(var_pair[1])
    ni, nj = PAR_NAMES[i], PAR_NAMES[j]
    fig, axes = plt.subplots(1, len(indices), figsize=(4.5 * len(indices), 4.5))
    circle = np.stack([np.cos(np.linspace(0, 2 * np.pi, 100)), np.sin(np.linspace(0, 2 * np.pi, 100))])
    for ax, idx in zip(np.atleast_1d(axes), indices):
        ax.scatter(post[var_pair[0]][idx].ravel(), post[var_pair[1]][idx].ravel(),
                  s=6, alpha=0.25, label="posterior draws")
        J_raw, theta0 = local_jacobian(workflow, val, idx, x_key, cond_keys=cond_keys)
        cov_pair = (J_raw @ J_raw.T)[np.ix_([i, j], [i, j])]
        eigvals, eigvecs = np.linalg.eigh(cov_pair)
        S_pair = np.sqrt(np.clip(eigvals, 0, None))
        center = theta0[[i, j]]
        ellipse = center[:, None] + (eigvecs * S_pair) @ circle
        ax.plot(ellipse[0], ellipse[1], color="crimson", linewidth=2, label=r"local Jacobian (1$\sigma$)")
        ax.scatter(*center, color="crimson", marker="x", zorder=5)
        ax.scatter(val[var_pair[0]][idx], val[var_pair[1]][idx], color="black", marker="*",
                  s=80, zorder=5, label="ground truth")
        ax.set_xlabel(ni); ax.set_ylabel(nj)
        ax.set_title(f"instance #{idx}")
        if bounds is not None:
            ax.set_xlim(bounds[var_pair[0]]); ax.set_ylim(bounds[var_pair[1]])
    np.atleast_1d(axes)[0].legend(fontsize=8, loc="best")
    fig.suptitle(f"Local identifiability (Jacobian SVD) -- {label}: {ni} vs {nj}", y=1.02)
    fig.tight_layout()
    return fig
