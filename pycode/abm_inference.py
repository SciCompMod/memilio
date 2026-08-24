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

SIM_DAYS = 30                    
REPORTING_DELAY_DAYS = 1         
OBSERVED_DAYS = SIM_DAYS - REPORTING_DELAY_DAYS  
SURVEY, DIAGNOSTIC = 0, 1        
N_CT_BINS = 41       


SOURCES_SURVEY = frozenset({SURVEY})
SOURCES_DIAGNOSTIC = frozenset({DIAGNOSTIC})
SOURCES_POOLED = frozenset({SURVEY, DIAGNOSTIC})

PARAM_KEYS = ["beta", "t_exposed", "time_presymptomatic", "time_asymptomatic_recovery", "symptom_prob",
             "quarantine_compliance", "care_seeking_prob"]
PAR_NAMES = [r"$\beta$", r"$t_\mathrm{exp}$", r"$t_\mathrm{presym}$", r"$t_\mathrm{asym}$", r"$p_\mathrm{sym}$",
            r"$p_\mathrm{quar}$", r"$p_\mathrm{care}$"]


PRIOR_BETA_S,   PRIOR_BETA_SCALE  = 0.5,  1.0   # median 1.0
PRIOR_TEXP_S,   PRIOR_TEXP_SCALE  = 0.4,  5.0   # median 5.0  (~95% in [2.3, 11])
PRIOR_TPRE_S,   PRIOR_TPRE_SCALE  = 0.5,  1.5   # median 1.5  (~95% in [0.56, 4.0])
PRIOR_TASYM_S,  PRIOR_TASYM_SCALE = 0.25, 8.0   # median 8.0  (~95% in [4.9, 13])
PRIOR_PSYM_MU,  PRIOR_PSYM_SIGMA  = np.log(0.6 / 0.4), 0.5  # logit-normal, median 0.6

PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA = np.log(0.6 / 0.4), 0.5  # logit-normal, median 0.6
PRIOR_CARE_MU,  PRIOR_CARE_SIGMA  = np.log(0.19 / 0.81), 0.5  # logit-normal, median 0.19

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
    z_care = rng.normal(PRIOR_CARE_MU, PRIOR_CARE_SIGMA)
    return {
        "beta":                       rng.lognormal(np.log(PRIOR_BETA_SCALE),  PRIOR_BETA_S),
        "t_exposed":                  rng.lognormal(np.log(PRIOR_TEXP_SCALE),  PRIOR_TEXP_S),
        "time_presymptomatic":        rng.lognormal(np.log(PRIOR_TPRE_SCALE),  PRIOR_TPRE_S),
        "time_asymptomatic_recovery": rng.lognormal(np.log(PRIOR_TASYM_SCALE), PRIOR_TASYM_S),
        "symptom_prob":               1.0 / (1.0 + np.exp(-z_sym)),  # expit(N(PSYM_MU, PSYM_SIGMA))
        "quarantine_compliance":      1.0 / (1.0 + np.exp(-z_qcomp)),  # expit(N(QCOMP_MU, QCOMP_SIGMA))
        "care_seeking_prob":          1.0 / (1.0 + np.exp(-z_care)),  # expit(N(CARE_MU, CARE_SIGMA))
    }


def _logitnorm_logpdf(p, mu, sigma):
    """log-density of a logit-normal at p in (0,1); -inf outside."""
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, -np.inf)
    m = (p > 0.0) & (p < 1.0)
    z = np.log(p[m] / (1.0 - p[m]))
    out[m] = stats.norm.logpdf(z, loc=mu, scale=sigma) - np.log(p[m]) - np.log1p(-p[m])
    return out


def _prior_logpdf_term(key, value):
    """log-density of ONE PARAM_KEYS component under prior_theta()'s independent prior -- the
    per-key building block prior_log_prob() sums over all PARAM_KEYS. Exposed standalone so
    estimate_mutual_information() can sum just a subset of terms when marginalizing (see
    `keep_keys` there)."""
    if key == "beta":
        return stats.lognorm.logpdf(value, s=PRIOR_BETA_S, scale=PRIOR_BETA_SCALE)
    elif key == "t_exposed":
        return stats.lognorm.logpdf(value, s=PRIOR_TEXP_S, scale=PRIOR_TEXP_SCALE)
    elif key == "time_presymptomatic":
        return stats.lognorm.logpdf(value, s=PRIOR_TPRE_S, scale=PRIOR_TPRE_SCALE)
    elif key == "time_asymptomatic_recovery":
        return stats.lognorm.logpdf(value, s=PRIOR_TASYM_S, scale=PRIOR_TASYM_SCALE)
    elif key == "symptom_prob":
        return _logitnorm_logpdf(value, PRIOR_PSYM_MU, PRIOR_PSYM_SIGMA)
    elif key == "quarantine_compliance":
        return _logitnorm_logpdf(value, PRIOR_QCOMP_MU, PRIOR_QCOMP_SIGMA)
    elif key == "care_seeking_prob":
        return _logitnorm_logpdf(value, PRIOR_CARE_MU, PRIOR_CARE_SIGMA)
    else:
        raise KeyError(f"no prior registered for {key!r}")


def prior_log_prob(**params):
    return sum(_prior_logpdf_term(k, params[k]) for k in PARAM_KEYS)



EPIDEMIC_PARAM_KEYS = [k for k in PARAM_KEYS if k not in ("quarantine_compliance", "care_seeking_prob")]


def prior_mean(n_mc=200_000, rng=None):
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
    z = rng.normal(PRIOR_CARE_MU, PRIOR_CARE_SIGMA, n_mc)
    means["care_seeking_prob"] = float(np.mean(1.0 / (1.0 + np.exp(-z))))
    return means


def n_events(period, sim_days=SIM_DAYS):
    return (sim_days - 1) // period + 1


def make_design(total_budget, period):
    event_budget_fraction = total_budget / n_events(period)
    return TestingBudget(event_budget_fraction=float(event_budget_fraction), test_period_days=int(period))


def total_budget_for_rate(rate_per_100k_per_day, sim_days=SIM_DAYS):

    return rate_per_100k_per_day / 100_000 * sim_days


def rate_per_100k_per_day(total_budget, sim_days=SIM_DAYS):
    return total_budget / sim_days * 100_000


def _assemble(draws, totals, periods, outs, sources=SOURCES_SURVEY, scale_by_tests=False):
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


def _marginal_kl_gaussian_plugin(samples, keep_keys, b, num_posterior_samples):
    diff = np.full((b, num_posterior_samples), np.nan)
    for i in range(b):
        theta_i = {k: samples[k].reshape(b, num_posterior_samples)[i] for k in keep_keys}
        u = np.stack([np.log1p(theta_i[k]) for k in keep_keys], axis=1)  # (num_posterior_samples, d)
        valid = np.all(np.isfinite(u), axis=1)
        if valid.sum() < len(keep_keys) + 2:  # not enough points for a full-rank covariance
            continue
        mean = u[valid].mean(axis=0)
        cov = np.cov(u[valid], rowvar=False)
        try:
            log_q_u = stats.multivariate_normal.logpdf(u, mean=mean, cov=cov, allow_singular=True)
        except Exception:
            continue
        log_p_theta = sum(_prior_logpdf_term(k, theta_i[k]) for k in keep_keys)
        log_p_u = log_p_theta + u.sum(axis=1)  # + log|Jacobian| of theta = expm1(u), elementwise
        row = log_q_u - log_p_u

        diff[i] = np.where(np.isfinite(row), row, np.nan)
    return diff


def estimate_mutual_information(workflow, val, x_key, num_posterior_samples=NUM_POSTERIOR_SAMPLES,
                                val_batch=32, cond_keys=("log_budget", "frequency"), keep_keys=None):
    """I(theta;X|d) for a validation set at a fixed design; the design (log_budget, frequency,
    carried in `val`) is passed as a condition to sample() and log_prob(). Pass cond_keys=()
    for a non-design-conditioned workflow.
    """
    keep_keys = PARAM_KEYS if keep_keys is None else list(keep_keys)
    marginal = keep_keys != PARAM_KEYS
    n_val = len(val["beta"])
    kls = np.empty(n_val)
    n_invalid, n_total = 0, 0
    for start in range(0, n_val, val_batch):
        idx = slice(start, min(start + val_batch, n_val))
        cond = {x_key: val[x_key][idx], **{k: val[k][idx] for k in cond_keys}}
        post = workflow.sample(conditions=cond, num_samples=num_posterior_samples)
        samples = {k: post[k].reshape(-1, 1) for k in PARAM_KEYS}
        b = val["beta"][idx].shape[0]
        if marginal:
            diff = _marginal_kl_gaussian_plugin(samples, keep_keys, b, num_posterior_samples)
        else:
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

    J_raw, theta_at_z = local_jacobian(workflow, val, idx, x_key, cond_keys, z)
    J_rel = np.diag(1.0 / theta_at_z) @ J_raw
    U, S, _ = np.linalg.svd(J_rel)
    return S, U, theta_at_z


def shared_bounds(val, posts, keys=PARAM_KEYS, pad=0.05, pct=(1, 99)):
    bounds = {}
    for k in keys:
        pooled = np.concatenate([val[k].ravel()] + [post[k].ravel() for post in posts.values()])
        lo, hi = np.percentile(pooled, pct)
        margin = pad * (hi - lo) if hi > lo else 1.0
        bounds[k] = (lo - margin, hi + margin)
    return bounds


def apply_recovery_bounds(fig, bounds, keys=PARAM_KEYS):

    for ax, k in zip(fig.axes, keys):
        ax.set_xlim(bounds[k]); ax.set_ylim(bounds[k])


def apply_pairgrid_bounds(g, bounds, keys=PARAM_KEYS):

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


def _kde_1d(samples, lo, hi, n_grid=200):
    """1-D KDE of `samples` evaluated on n_grid points over [lo, hi]; None if the fit fails
    (e.g. near-zero variance)."""
    try:
        kde = stats.gaussian_kde(samples)
    except Exception:
        return None
    xs = np.linspace(lo, hi, n_grid)
    return xs, kde(xs)


def _kde_2d(x, y, xlo, xhi, ylo, yhi, n_grid=60):
    """2-D KDE of (x, y) evaluated on an n_grid x n_grid mesh over [xlo,xhi] x [ylo,yhi]; None
    if the fit fails (e.g. a near-degenerate/collapsed posterior makes the sample covariance
    singular)."""
    try:
        kde = stats.gaussian_kde(np.vstack([x, y]))
    except Exception:
        return None
    xx, yy = np.mgrid[xlo:xhi:complex(n_grid), ylo:yhi:complex(n_grid)]
    zz = kde(np.vstack([xx.ravel(), yy.ravel()])).reshape(xx.shape)
    return xx, yy, zz


def plot_posterior_comparison_grid(posts, val, idx, keys=PARAM_KEYS, names=PAR_NAMES,
                                   labels=("ct", "bin"), colors=("tab:blue", "tab:orange"),
                                   bounds=None, n_levels=4):
    n = len(keys)
    fig, axes = plt.subplots(n, n, figsize=(2.3 * n, 2.3 * n))
    for row in range(n):
        for col in range(n):
            ax = axes[row, col]
            if col > row:
                ax.axis("off")
                continue
            kr, kc = keys[row], keys[col]
            if row == col:
                for label, color in zip(labels, colors):
                    s = posts[label][kr][idx].ravel()
                    lo, hi = bounds[kr] if bounds is not None else (s.min(), s.max())
                    fit = _kde_1d(s, lo, hi)
                    if fit is not None:
                        ax.plot(*fit, color=color, label=label)
                ax.axvline(val[kr][idx], color="black", linestyle="--", linewidth=1)
            else:
                for label, color in zip(labels, colors):
                    x = posts[label][kc][idx].ravel()
                    y = posts[label][kr][idx].ravel()
                    xlo, xhi = bounds[kc] if bounds is not None else (x.min(), x.max())
                    ylo, yhi = bounds[kr] if bounds is not None else (y.min(), y.max())
                    fit = _kde_2d(x, y, xlo, xhi, ylo, yhi)
                    if fit is not None:
                        ax.contour(*fit, levels=n_levels, colors=color, linewidths=1)
                    else:
                        ax.scatter(x, y, s=4, alpha=0.2, color=color)
                ax.scatter(val[kc][idx], val[kr][idx], color="black", marker="*", s=70, zorder=5)
                if bounds is not None:
                    ax.set_xlim(bounds[kc]); ax.set_ylim(bounds[kr])
            if row == n - 1:
                ax.set_xlabel(names[col])
            else:
                ax.set_xticklabels([])
            if col == 0:
                ax.set_ylabel(names[row])
            else:
                ax.set_yticklabels([])
    handles = [plt.Line2D([0], [0], color=c, label=l) for l, c in zip(labels, colors)]
    handles.append(plt.Line2D([0], [0], color="black", marker="*", linestyle="none", label="ground truth"))
    fig.legend(handles=handles, loc="upper right", fontsize=9)
    fig.suptitle(f"Posterior comparison ({' vs '.join(labels)}) -- instance #{idx}", y=1.0)
    return fig
