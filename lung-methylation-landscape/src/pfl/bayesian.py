"""Conjugate Bayesian regression with effective-sample-size power weights."""
from __future__ import annotations

import numpy as np
from scipy.stats import invgamma
from scipy.stats import t as student_t


def effective_weights(weights):
    """Scale nonnegative terminal weights to Kish effective sample size."""
    weights = np.asarray(weights, float)
    weights = np.clip(weights, 0, None)
    if np.all(weights == 0):
        raise ValueError("All terminal conditioning weights are zero")
    kish = float(weights.sum() ** 2 / (np.sum(weights * weights) + 1e-15))
    scaled = weights * (kish / (weights.sum() + 1e-15))
    return scaled, kish


def nig_fit(X, y, w=None, *, n_draw, prior_variance, a0, b0, rng):
    """Fit a normal-inverse-gamma model and sample posterior coefficients."""
    X = np.asarray(X, float)
    y = np.asarray(y, float)
    valid = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
    X, y = X[valid], y[valid]

    if w is None:
        weights, kish = np.ones(len(y)), float(len(y))
    else:
        weights, kish = effective_weights(np.asarray(w, float)[valid])

    square_root = np.sqrt(weights)
    Xw, yw = X * square_root[:, None], y * square_root
    n_eff, n_coefficients = float(weights.sum()), X.shape[1]

    prior_covariance = np.eye(n_coefficients) * prior_variance
    prior_covariance[0, 0] = 100.0
    prior_precision = np.linalg.inv(prior_covariance)
    posterior_precision = prior_precision + Xw.T @ Xw
    covariance = np.linalg.inv(posterior_precision)
    mean = covariance @ (Xw.T @ yw)

    an = a0 + n_eff / 2.0
    bn = b0 + 0.5 * (float(yw @ yw) - float(mean @ (posterior_precision @ mean)))
    bn = max(float(bn), 1e-12)
    df = 2 * an
    scale = np.sqrt((bn / an) * np.diag(covariance))
    lower = student_t.ppf(0.025, df, loc=mean, scale=scale)
    upper = student_t.ppf(0.975, df, loc=mean, scale=scale)
    p_positive = 1 - student_t.cdf(0, df, loc=mean, scale=scale)
    p_negative = student_t.cdf(0, df, loc=mean, scale=scale)

    sigma_squared = invgamma.rvs(an, scale=bn, size=n_draw, random_state=rng)
    L = np.linalg.cholesky(covariance)
    standard_normal = rng.normal(size=(n_draw, n_coefficients))
    draws = mean[None, :] + (standard_normal @ L.T) * np.sqrt(sigma_squared)[:, None]
    return {
        'n': len(y), 'effective_n': kish, 'p': n_coefficients,
        'mn': mean, 'Vn': covariance, 'an': an, 'bn': bn, 'df': df,
        'scale': scale, 'qlo': lower, 'qhi': upper,
        'pgt': p_positive, 'plt': p_negative, 'draws': draws,
    }
