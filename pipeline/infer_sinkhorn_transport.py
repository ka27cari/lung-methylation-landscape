#!/usr/bin/env python3
"""Add an entropy-regularized transport route layer to the frozen landscape.

This script deliberately does *not* redefine terminal states or the primary
landscape.  Terminal identities remain those obtained from the exact partial
OT sensitivity rule in pipeline/fit_microstate_landscape.py.  Here, balanced entropic
Sinkhorn OT is used only as a smoother downstream representation of transport
between the already-frozen source and prediagnostic state occupancies.

Primary route representation
----------------------------
* balanced Sinkhorn OT, epsilon = config/model.yaml:sinkhorn_route.epsilon
* antisymmetric net current J_ij = max(G_ij - G_ji, 0)
* visualization support = smallest descending set of positive J edges carrying
  the configured fraction (default 95%) of total net-current mass

Outputs are written to cache/model and analysis/tables.  A small epsilon
sensitivity audit is also emitted, but it does not participate in model
selection.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import yaml
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[1] / "rebuild"
MOD = ROOT / 'cache/model'
TAB = ROOT / 'analysis/tables'
MOD.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

CFG = yaml.safe_load((ROOT.parent / 'config/model.yaml').read_text())
RCFG = CFG.get('sinkhorn_route', {})
EPS = float(RCFG.get('epsilon', 0.10))
EPS_SENS = [float(x) for x in RCFG.get('epsilon_sensitivity', [0.05, 0.10, 0.20, 0.40])]
RETAIN = float(RCFG.get('net_current_mass_fraction', 0.95))
MAX_ITER = int(RCFG.get('max_iter', 50000))
TOL = float(RCFG.get('tol', 1e-13))


def sinkhorn(a, b, C, eps, n_iter=50000, tol=1e-8):
    """Balanced entropy-regularized OT using classic matrix scaling.

    C is already median-normalized, matching the exact partial-OT cost scale.
    For the present 16-state problem this direct implementation is stable and
    keeps the build independent of an additional POT dependency.  Convergence
    is assessed on the transported row/column marginals, which is the relevant
    numerical criterion for the coupling used downstream.
    """
    K = np.exp(-C / eps)
    u = np.ones_like(a, dtype=float)
    v = np.ones_like(b, dtype=float)
    tiny = 1e-300
    converged = False
    it = 0
    marginal_err = np.inf
    for it in range(1, n_iter + 1):
        u = np.divide(a, K @ v + tiny)
        v = np.divide(b, K.T @ u + tiny)
        if it == 1 or it % 25 == 0:
            G = (u[:, None] * K) * v[None, :]
            row_err = float(np.max(np.abs(G.sum(axis=1) - a)))
            col_err = float(np.max(np.abs(G.sum(axis=0) - b)))
            marginal_err = max(row_err, col_err)
            if marginal_err < tol:
                converged = True
                break
    G = (u[:, None] * K) * v[None, :]
    row_err = float(np.max(np.abs(G.sum(axis=1) - a)))
    col_err = float(np.max(np.abs(G.sum(axis=0) - b)))
    marginal_err = max(row_err, col_err)
    if marginal_err < tol:
        converged = True
    return G, dict(iterations=it, converged=converged, marginal_error=marginal_err,
                   max_row_marginal_error=row_err, max_col_marginal_error=col_err)


def net_current(G):
    J = np.maximum(G - G.T, 0.0)
    np.fill_diagonal(J, 0.0)
    return J


def retain_top_mass(J, fraction=0.95):
    positive = [(float(J[i, j]), i, j)
                for i in range(J.shape[0])
                for j in range(J.shape[1]) if J[i, j] > 0]
    positive.sort(reverse=True)
    total = sum(x[0] for x in positive)
    out = np.zeros_like(J)
    cumulative = 0.0
    if total <= 0:
        return out, 0, 0.0
    for mass, i, j in positive:
        out[i, j] = mass
        cumulative += mass
        if cumulative / total >= fraction:
            break
    return out, int(np.sum(out > 0)), cumulative / total


def absorption_from_plan(G, terms):
    """Terminal absorption affinities on a directed route matrix.

    Diagonal mass is removed, each nonterminal row is normalized to transition
    probabilities, and the frozen LP-defined terminals are made absorbing.
    Rows with no terminal-reachable outgoing route retain zero affinity.
    """
    n = G.shape[0]
    terms = list(map(int, terms))
    termset = set(terms)
    non = [s for s in range(n) if s not in termset]
    ni = {s: i for i, s in enumerate(non)}
    ti = {s: i for i, s in enumerate(terms)}
    Q = np.zeros((len(non), len(non)), dtype=float)
    R = np.zeros((len(non), len(terms)), dtype=float)
    for s in non:
        row = np.asarray(G[s], float).copy()
        row[s] = 0.0
        tot = float(row.sum())
        if tot <= 1e-15:
            continue
        row /= tot
        for j, val in enumerate(row):
            if val <= 0:
                continue
            if j in termset:
                R[ni[s], ti[j]] += val
            elif j in ni:
                Q[ni[s], ni[j]] += val
    I = np.eye(len(non))
    try:
        N = np.linalg.inv(I - Q)
    except np.linalg.LinAlgError:
        N = np.linalg.pinv(I - Q)
    B = N @ R
    F = N @ B
    H = np.divide(F, B, out=np.full_like(F, np.nan), where=B > 1e-12)
    A = np.zeros((n, len(terms)), dtype=float)
    HH = np.full((n, len(terms)), np.nan, dtype=float)
    for s in non:
        A[s] = B[ni[s]]
        HH[s] = H[ni[s]]
    for t in terms:
        A[t, ti[t]] = 1.0
        HH[t, ti[t]] = 0.0
    return A, HH


roles = pd.read_csv(TAB / 'FINAL_STATE_ROLE_CLASSIFICATION.csv')
roles['state'] = roles['state'].astype(int)
cent = pd.read_csv(MOD / 'STATE_CENTROIDS.csv').sort_values('state')
state_order = cent['state'].astype(int).to_numpy()
if not np.array_equal(state_order, np.arange(len(state_order))):
    raise RuntimeError('STATE_CENTROIDS.csv must contain contiguous state IDs 0..K-1')
program_cols = [c for c in cent.columns if c.startswith('program_')]
C0 = cdist(cent[program_cols].to_numpy(float), cent[program_cols].to_numpy(float), 'sqeuclidean')
pos = C0[C0 > 0]
C = C0 / (np.median(pos) if len(pos) else 1.0)

roles = roles.sort_values('state').reset_index(drop=True)
a = roles['control_occupancy'].to_numpy(float)
b = roles['prediagnostic_occupancy'].to_numpy(float)
# Numerical normalization only; the generated occupancies should already sum to 1.
a = a / a.sum()
b = b / b.sum()
terms = sorted(roles.loc[roles['final_microstate_terminal'].astype(bool), 'state'].astype(int).tolist())
K = len(a)

G, diag = sinkhorn(a, b, C, EPS, MAX_ITER, TOL)
J = net_current(G)
J95, n95, retained_actual = retain_top_mass(J, RETAIN)
A_raw, H_raw = absorption_from_plan(G, terms)
A_net95, H_net95 = absorption_from_plan(J95, terms)

state_names = [f'S{i}' for i in range(K)]
pd.DataFrame(G, index=state_names, columns=state_names).to_csv(MOD / 'SINKHORN_OT_PLAN.csv')
pd.DataFrame(J, index=state_names, columns=state_names).to_csv(MOD / 'SINKHORN_NET_CURRENT.csv')
pd.DataFrame(J95, index=state_names, columns=state_names).to_csv(MOD / 'SINKHORN_NET_CURRENT_95.csv')

# Long edge table is convenient for R / plotting and keeps the retained flag explicit.
edge_rows = []
for i in range(K):
    for j in range(K):
        if i == j or J[i, j] <= 0:
            continue
        edge_rows.append(dict(
            source_state=i,
            target_state=j,
            sinkhorn_mass=float(G[i, j]),
            reverse_sinkhorn_mass=float(G[j, i]),
            net_current=float(J[i, j]),
            retained_net95=bool(J95[i, j] > 0),
        ))
edge_df = pd.DataFrame(edge_rows).sort_values('net_current', ascending=False)
edge_df.to_csv(TAB / 'SINKHORN_NET_CURRENT_EDGES.csv', index=False)

# Route affinity tables use the frozen LP terminals; no terminal reclassification occurs here.
route_rows = []
for s in range(K):
    for jj, t in enumerate(terms):
        route_rows.append(dict(
            state=s,
            terminal_state=t,
            raw_sinkhorn_absorption_affinity=float(A_raw[s, jj]),
            raw_sinkhorn_expected_steps=float(H_raw[s, jj]) if np.isfinite(H_raw[s, jj]) else np.nan,
            net95_absorption_affinity=float(A_net95[s, jj]),
            net95_expected_steps=float(H_net95[s, jj]) if np.isfinite(H_net95[s, jj]) else np.nan,
        ))
pd.DataFrame(route_rows).to_csv(TAB / 'SINKHORN_TERMINAL_ROUTE_AFFINITY.csv', index=False)

# Epsilon sensitivity is descriptive only: edge count, marginal error, and preferred terminal stability.
sens_rows = []
for eps in EPS_SENS:
    Ge, de = sinkhorn(a, b, C, eps, MAX_ITER, TOL)
    Je = net_current(Ge)
    Jkeep, nkeep, frac = retain_top_mass(Je, RETAIN)
    Ae, _ = absorption_from_plan(Ge, terms)
    for s in range(K):
        if s in terms:
            pref = s
            pref_aff = 1.0
        elif Ae.shape[1] and float(Ae[s].sum()) > 1e-12:
            jj = int(np.argmax(Ae[s]))
            pref = int(terms[jj])
            pref_aff = float(Ae[s, jj])
        else:
            pref = np.nan
            pref_aff = np.nan
        sens_rows.append(dict(
            epsilon=float(eps),
            state=int(s),
            role=str(roles.loc[roles.state.eq(s), 'macrostate_role'].iloc[0]),
            preferred_terminal=pref,
            preferred_affinity=pref_aff,
            raw_positive_offdiag_edges=int(np.sum((Ge > 0) & (~np.eye(K, dtype=bool)))),
            positive_net_current_edges=int(np.sum(Je > 0)),
            retained_net_current_edges=int(nkeep),
            retained_net_current_fraction=float(frac),
            max_row_marginal_error=float(de['max_row_marginal_error']),
            max_col_marginal_error=float(de['max_col_marginal_error']),
        ))
pd.DataFrame(sens_rows).to_csv(TAB / 'SINKHORN_ROUTE_SENSITIVITY.csv', index=False)

summary = dict(
    purpose='downstream route representation only; terminal identities remain frozen from source occupancy-gain + primary exact partial-LP95 discovery',
    epsilon=EPS,
    epsilon_sensitivity=EPS_SENS,
    cost='squared Euclidean state-centroid distance normalized by the median positive pairwise cost',
    balanced_marginals='control and prediagnostic state occupancies, each normalized to total mass 1',
    net_current='J_ij=max(G_ij-G_ji,0), diagonal set to 0',
    net_current_mass_fraction_requested=RETAIN,
    net_current_mass_fraction_retained=float(retained_actual),
    raw_positive_offdiag_edges=int(np.sum((G > 0) & (~np.eye(K, dtype=bool)))),
    positive_net_current_edges=int(np.sum(J > 0)),
    retained_net_current_edges=int(n95),
    terminal_states=terms,
    sinkhorn_iterations=int(diag['iterations']),
    sinkhorn_converged=bool(diag['converged']),
    max_row_marginal_error=float(diag['max_row_marginal_error']),
    max_col_marginal_error=float(diag['max_col_marginal_error']),
)
json.dump(summary, open(TAB / 'SINKHORN_ROUTE_LAYER_SPEC.json', 'w'), indent=2)

# Append, rather than overwrite, route-layer metadata in the canonical model spec.
spec_path = TAB / 'FINAL_MODEL_SPEC.json'
if spec_path.exists():
    spec = json.load(open(spec_path))
    spec['route_transport'] = dict(
        method='balanced_entropic_sinkhorn_net_current',
        role='downstream route representation; not used to define terminals',
        epsilon=EPS,
        net_current_mass_fraction=RETAIN,
        terminal_definition='frozen terminals from source occupancy-gain + primary exact partial-LP95 destination rule; LP90/99 are sensitivity only',
    )
    json.dump(spec, open(spec_path, 'w'), indent=2)

print(json.dumps(summary, indent=2))
