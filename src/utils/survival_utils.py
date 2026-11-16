# =============================================================================
# src/utils/survival_utils.py
# Shared survival utility functions.
# =============================================================================

import numpy as np
import pandas as pd


def determine_surv_prob(surv, t):
    """Interpolate survival probability at time t from a survival DataFrame."""
    for i, t2 in enumerate(surv.index):
        if t2 > t:
            if i == 0:
                result = (t / t2) * surv.iloc[i] + (t2 - t) / t2 * 1
            else:
                t1 = surv.index.values[i - 1]
                result = ((t - t1) / (t2 - t1) * surv.iloc[i] +
                          (t2 - t) / (t2 - t1) * surv.iloc[i - 1])
            break
        else:
            result = surv.iloc[i]
    return result


def surv_at_horizons(surv, times_to_predict):
    """
    Evaluate survival probabilities at given time horizons.
    surv: DataFrame with index = time points, columns = subjects.
    Returns DataFrame with columns = ["P25", "P5", "P75"].
    """
    pred = pd.DataFrame(
        [list(determine_surv_prob(surv, t)) for t in times_to_predict]
    ).T
    pred.columns = ["P25", "P5", "P75"]
    return pred


def rmst_from_surv(surv, horizon):
    """
    Compute RMST up to horizon for each subject via trapezoidal rule.
    surv: DataFrame with index = time points, columns = subjects.
    Returns a Series of RMST values per subject.
    """
    times = np.array(surv.index, dtype=float)
    mask  = times <= horizon
    t_use = np.concatenate([[0], times[mask], [horizon]])

    S = surv.values.T  # shape: (n_subjects, n_times)

    # interpolate S at horizon
    if mask.sum() == 0:
        s_at_horizon = np.ones(S.shape[0])
    else:
        s_at_horizon = S[:, mask][:, -1]

    S_use = np.concatenate([
        np.ones((S.shape[0], 1)),
        S[:, mask],
        s_at_horizon[:, None]
    ], axis=1)

    dt   = np.diff(t_use)
    avgS = (S_use[:, :-1] + S_use[:, 1:]) / 2
    return pd.Series(avgS @ dt, index=surv.columns)