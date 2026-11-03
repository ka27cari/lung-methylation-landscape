"""Control-based transfer of frozen methylation landscape coordinates."""
from __future__ import annotations

import numpy as np


def align_controls(external, source_controls, external_controls):
    """Match external control medians/IQRs to the source control reference."""
    external = np.asarray(external, dtype=float)
    source_controls = np.asarray(source_controls, dtype=float)
    external_controls = np.asarray(external_controls, dtype=float)
    if any(x.ndim != 2 for x in (external, source_controls, external_controls)):
        raise ValueError("Expected 2D program score matrices")
    if len({x.shape[1] for x in (external, source_controls, external_controls)}) != 1:
        raise ValueError("Program dimensions do not match")
    if not len(source_controls) or not len(external_controls):
        raise ValueError("Both reference control groups must be nonempty")

    ext_iqr = np.percentile(external_controls, 75, axis=0) - np.percentile(
        external_controls, 25, axis=0
    )
    src_iqr = np.percentile(source_controls, 75, axis=0) - np.percentile(
        source_controls, 25, axis=0
    )
    if (ext_iqr <= 0).any() or not np.isfinite(ext_iqr).all():
        raise ValueError("External reference contains zero or invalid IQR")
    return (
        (external - np.median(external_controls, axis=0))
        * (src_iqr / ext_iqr)
        + np.median(source_controls, axis=0)
    )


def state_posterior(scores, centroids, state_scales, temperature):
    """Compute frozen soft state memberships for normalized program scores."""
    x = np.asarray(scores, dtype=float)
    c = np.asarray(centroids, dtype=float)
    scales = np.asarray(state_scales, dtype=float)
    if x.ndim != 2 or c.ndim != 2 or x.shape[1] != c.shape[1]:
        raise ValueError("Invalid sample/centroid dimensions")
    if len(scales) != len(c) or not np.all(scales > 0) or temperature <= 0:
        raise ValueError("Invalid posterior scale or temperature")
    distance_squared = np.sum((x[:, None, :] - c[None, :, :]) ** 2, axis=2)
    logits = -distance_squared / (2.0 * float(temperature) * scales[None, :])
    logits -= logits.max(axis=1, keepdims=True)
    weights = np.exp(logits)
    return weights / weights.sum(axis=1, keepdims=True)


def terminal_mass(scores, centroids, state_scales, temperature, terminal_states):
    """Sum frozen posterior membership over the selected terminal states."""
    posterior = state_posterior(scores, centroids, state_scales, temperature)
    terminals = _terminal_indices(terminal_states, posterior.shape[1])
    return posterior[:, terminals].sum(axis=1)


def _terminal_indices(terminal_states, n_states):
    terminals = np.asarray(terminal_states, dtype=int)
    if not len(terminals) or (terminals < 0).any() or (terminals >= n_states).any():
        raise ValueError("Invalid terminal-state indices")
    return terminals


def frozen_projection(external_scores, source_scores, source_control_mask,
                      external_control_mask, centroids, state_scales,
                      temperature, terminal_states):
    """Return calibrated/raw program scores and their frozen posteriors."""
    source = np.asarray(source_scores, dtype=float)
    external = np.asarray(external_scores, dtype=float)
    src_ctrl = np.asarray(source_control_mask, dtype=bool)
    ext_ctrl = np.asarray(external_control_mask, dtype=bool)
    if len(source) != len(src_ctrl) or len(external) != len(ext_ctrl):
        raise ValueError("Control masks must match their sample matrices")

    aligned = align_controls(external, source[src_ctrl], external[ext_ctrl])
    posterior_raw = state_posterior(external, centroids, state_scales, temperature)
    posterior = state_posterior(aligned, centroids, state_scales, temperature)
    terminals = _terminal_indices(terminal_states, posterior.shape[1])
    return {
        "raw_scores": external,
        "scores": aligned,
        "raw_posterior": posterior_raw,
        "posterior": posterior,
        "raw_terminal_mass": posterior_raw[:, terminals].sum(axis=1),
        "terminal_mass": posterior[:, terminals].sum(axis=1),
    }
