"""Two-sample Mendelian randomisation estimators on harmonised summary statistics.

Conventions follow TwoSampleMR so results are comparable with the R package:
Wald ratio SE is first order (se_out / |b_exp|); IVW fixed effect uses
inverse-variance weights 1/se_ratio^2; the multiplicative random-effects IVW
inflates the fixed-effect SE by sqrt(max(1, Q / (k - 1))).
"""
from __future__ import annotations

import numpy as np
from scipy import stats


def _res(method, beta, se, k, **extra):
    beta = float(beta)
    se = float(se)
    p = float(2 * stats.norm.sf(abs(beta / se))) if se > 0 else np.nan
    return {"method": method, "n_snp": int(k), "beta": beta, "se": se, "p": p,
            "OR": float(np.exp(beta)), "OR_lci": float(np.exp(beta - 1.959964 * se)),
            "OR_uci": float(np.exp(beta + 1.959964 * se)), **extra}


def wald_ratio(b_exp, se_exp, b_out, se_out, second_order=False):
    b_exp, se_exp, b_out, se_out = map(float, (b_exp, se_exp, b_out, se_out))
    beta = b_out / b_exp
    se = abs(se_out / b_exp)
    if second_order:
        se = np.sqrt(se_out ** 2 / b_exp ** 2 + b_out ** 2 * se_exp ** 2 / b_exp ** 4)
    return _res("Wald ratio", beta, se, 1)


def ivw(b_exp, se_exp, b_out, se_out):
    """Fixed-effect and multiplicative random-effects IVW with Cochran's Q."""
    b_exp, se_exp, b_out, se_out = (np.asarray(x, float) for x in (b_exp, se_exp, b_out, se_out))
    k = len(b_exp)
    ratio = b_out / b_exp
    se_ratio = np.abs(se_out / b_exp)
    w = 1.0 / se_ratio ** 2
    beta = np.sum(w * ratio) / np.sum(w)
    se_fe = np.sqrt(1.0 / np.sum(w))
    q = float(np.sum(w * (ratio - beta) ** 2))
    df = k - 1
    q_p = float(stats.chi2.sf(q, df)) if df > 0 else np.nan
    phi = q / df if df > 0 else np.nan
    fe = _res("IVW (fixed effect)", beta, se_fe, k, Q=q, Q_df=df, Q_p=q_p)
    se_re = se_fe * np.sqrt(max(1.0, phi)) if df > 0 else se_fe
    re = _res("IVW (multiplicative random effects)", beta, se_re, k, Q=q, Q_df=df, Q_p=q_p, phi=phi)
    return fe, re


def mr_egger(b_exp, se_exp, b_out, se_out):
    """MR-Egger regression, orienting variants so exposure effects are positive."""
    b_exp, se_exp, b_out, se_out = (np.asarray(x, float) for x in (b_exp, se_exp, b_out, se_out))
    k = len(b_exp)
    if k < 3:
        return None
    sign = np.sign(b_exp)
    x = np.abs(b_exp)
    y = b_out * sign
    w = 1.0 / se_out ** 2
    X = np.column_stack([np.ones(k), x])
    W = np.diag(w)
    xtwx_inv = np.linalg.inv(X.T @ W @ X)
    coef = xtwx_inv @ X.T @ W @ y
    resid = y - X @ coef
    sigma2 = float(np.sum(w * resid ** 2) / (k - 2))
    cov = xtwx_inv * max(1.0, sigma2)  # residual SE floored at 1, as in TwoSampleMR
    se = np.sqrt(np.diag(cov))
    slope = _res("MR-Egger", coef[1], se[1], k)
    t = coef[1] / se[1]
    slope["p"] = float(2 * stats.t.sf(abs(t), k - 2))
    t0 = coef[0] / se[0]
    slope.update(intercept=float(coef[0]), intercept_se=float(se[0]),
                 intercept_p=float(2 * stats.t.sf(abs(t0), k - 2)))
    return slope


def weighted_median(b_exp, se_exp, b_out, se_out, n_boot=1000, seed=1):
    """Weighted median estimator (Bowden et al. 2016) with bootstrap SE."""
    b_exp, se_exp, b_out, se_out = (np.asarray(x, float) for x in (b_exp, se_exp, b_out, se_out))
    k = len(b_exp)
    if k < 3:
        return None

    def wm(ratio, w):
        order = np.argsort(ratio)
        r, ww = ratio[order], w[order] / np.sum(w)
        cum = np.cumsum(ww) - 0.5 * ww
        below = np.max(np.where(cum < 0.5)[0]) if np.any(cum < 0.5) else 0
        if below + 1 >= len(r):
            return r[-1]
        return r[below] + (r[below + 1] - r[below]) * (0.5 - cum[below]) / (cum[below + 1] - cum[below])

    ratio = b_out / b_exp
    w = (b_exp / se_out) ** 2
    est = wm(ratio, w)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        be = rng.normal(b_exp, se_exp)
        bo = rng.normal(b_out, se_out)
        boots.append(wm(bo / be, (be / se_out) ** 2))
    return _res("Weighted median", est, np.std(boots, ddof=1), k)


def run_mr(inst):
    """Run the appropriate estimators on a harmonised instrument table.

    ``inst`` needs columns beta_exp, se_exp, beta_out, se_out. Returns
    (primary, [all results]); the primary estimate is the Wald ratio for one
    instrument and the multiplicative random-effects IVW otherwise.
    """
    args = (inst.beta_exp.values, inst.se_exp.values, inst.beta_out.values, inst.se_out.values)
    if len(inst) == 1:
        w = wald_ratio(*(a[0] for a in args))
        return w, [w]
    fe, re = ivw(*args)
    out = [fe, re]
    for f in (mr_egger, weighted_median):
        r = f(*args)
        if r is not None:
            out.append(r)
    return re, out
