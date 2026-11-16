# =============================================================================
# scripts/run_nn.py
#
# Fit CoxCC or CoxTime (pycox feedforward neural networks) on simulated data
# from setting 1 or 2.  Edit MODEL_NAME and SETTING below, then run:
#
#   python scripts/run_nn.py
#
# from the repository root.  Hyperparameters are selected with Optuna (TPE).
#
# Output per replication j  (mirrors 02_run_models.R):
#   results/sim_setting<k>/<MODEL>/<scenario>/simdata_<j>/
#     best_params.csv        — selected hyperparameters
#     pred.csv               — full survival curve on test set (rows = time)
#     pred_horizons.csv      — factual survival probabilities at P25, P5, P75
#     rmst_horizons.csv      — factual RMST at P25, P5, P75
#     pred_alt.csv           — counterfactual survival curve
#     pred_horizons_alt.csv  — counterfactual survival probabilities
#     rmst_horizons_alt.csv  — counterfactual RMST
#
# Paper hyperparameter search space (Supplemental Material Table 1):
#   n_layers      : 1–4
#   neurons       : 32–256 (powers of 2)
#   dropout       : 0.0–0.5
#   l2            : 1e-5–1e-2 (log-uniform)
#   learning_rate : 1e-4–1e-2 (log-uniform)
#   batch_size    : 128–512 (powers of 2)
#   optimizer     : Adam, SGD+cosine-annealing, Adam-AMSGrad
#   activation    : ReLU, Tanh
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd
import optuna
optuna.logging.set_verbosity(optuna.logging.WARNING)

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "src"))

from models.deep_survival import make_mapper, objective, fit_model, predict_surv
from utils.survival_utils import surv_at_horizons, rmst_from_surv

# =============================================================================
# USER SETTINGS — edit here
# =============================================================================

MODEL_NAME = "CoxTime"   # "CoxTime" or "CoxCC"
SETTING    = 1           # 1 = Weibull/Biospear  |  2 = AFT/Friedman
N_TRIALS   = 200         # Optuna trials per replication (paper uses 200)
ITERS      = list(range(1, 101))   # replications 1 … 100

# Must match values used in 01_generate_data.R
S1_N = 2000
S2_N = 2000

S1_ALPHA = abs(-0.288)   # |log(0.75)|
S1_BETA  = 0.288
S1_GAMMA = 0.9           # |−0.9|

S2_A_MIN   = 0.2         # |−0.2|
S2_A_MAX   = 0.3
S2_NVARTOT = 5
S2_NVARACT = 5
S2_L       = 10

# =============================================================================
# Resolve scenario label and paths
# =============================================================================
if SETTING == 1:
    setting_root = os.path.join(REPO_ROOT, "data", "sim_setting1")
    sim_details  = os.listdir(setting_root)[0]
    data_root    = os.path.join(setting_root, sim_details)
    res_root     = os.path.join(REPO_ROOT, "results", "sim_setting1", MODEL_NAME)
elif SETTING == 2:
    setting_root = os.path.join(REPO_ROOT, "data", "sim_setting2")
    sim_details  = os.listdir(setting_root)[0]
    data_root    = os.path.join(setting_root, sim_details)
    res_root     = os.path.join(REPO_ROOT, "results", "sim_setting2", MODEL_NAME)
else:
    raise ValueError("SETTING must be 1 or 2.")

os.makedirs(res_root, exist_ok=True)
print(f"Setting {SETTING} | model: {MODEL_NAME} | scenario: {sim_details}")

# =============================================================================
# Config container
# =============================================================================
class Config:
    """Holds hyperparameters for the neural network."""
    pass

# =============================================================================
# Main loop
# =============================================================================
for j in ITERS:

    print(f"  replication {j} / {max(ITERS)}")

    dir_res = os.path.join(res_root, sim_details, f"simdata_{j}")
    os.makedirs(dir_res, exist_ok=True)

    data_dir = os.path.join(data_root, "all", f"simdata_{j}")
    df_train = pd.read_csv(os.path.join(data_dir, "train.csv"), index_col=0)
    df_test  = pd.read_csv(os.path.join(data_dir, "test.csv"),  index_col=0)

    times_path       = os.path.join(data_root, "surv_th", f"simdata_{j}", "times.csv")
    times            = pd.read_csv(times_path)
    times_to_predict = [times["P25"].item(), times["P5"].item(), times["P75"].item()]

    config      = Config()
    config.name = MODEL_NAME
    df_mapper   = make_mapper(df_train)

    db_path = os.path.join(dir_res, f"{MODEL_NAME}_simdata_{j}.db")
    study   = optuna.create_study(
        study_name     = f"{MODEL_NAME}_{sim_details}_simdata_{j}",
        storage        = f"sqlite:///{db_path}",
        sampler        = optuna.samplers.TPESampler(seed=j),
        direction      = "minimize",
        load_if_exists = True,
    )
    study.optimize(
        lambda trial: objective(trial, df_train, df_mapper, config),
        n_trials=N_TRIALS,
        show_progress_bar=False,
    )

    best = study.best_trial.params
    config.acti_func  = best["activation"]
    config.batch_size = best["batch_size"]
    config.layers     = best["n_layers"]
    config.lr         = best["learning_rate"]
    config.neurons    = best["neurons"]
    config.optim      = best["optimizer"]
    config.pen_l2     = best["l2"]
    config.dr         = best["dropout"]

    pd.DataFrame([best]).to_csv(os.path.join(dir_res, "best_params.csv"), index=False)
    df_val   = df_train.sample(frac=0.2, random_state=j)
    df_train_fit = df_train.drop(df_val.index)
    model, df_mapper = fit_model(config, df_train_fit, df_val, df_mapper)

    surv          = predict_surv(model, df_test, df_mapper)
    pred_horizons = surv_at_horizons(surv, times_to_predict)
    rmst_horizons = pd.DataFrame({
        col: rmst_from_surv(surv, t)
        for col, t in zip(["P25", "P5", "P75"], times_to_predict)
    })

    surv.T.to_csv(          os.path.join(dir_res, "pred.csv"),              index=False)
    pred_horizons.to_csv(   os.path.join(dir_res, "pred_horizons.csv"),     index=False)
    rmst_horizons.to_csv(   os.path.join(dir_res, "rmst_horizons.csv"),     index=False)

    df_test_alt        = df_test.copy()
    df_test_alt["Trt"] = -df_test_alt["Trt"]
    surv_alt           = predict_surv(model, df_test_alt, df_mapper)
    pred_alt_horizons  = surv_at_horizons(surv_alt, times_to_predict)
    rmst_alt_horizons  = pd.DataFrame({
        col: rmst_from_surv(surv_alt, t)
        for col, t in zip(["P25", "P5", "P75"], times_to_predict)
    })

    surv_alt.T.to_csv(        os.path.join(dir_res, "pred_alt.csv"),             index=False)
    pred_alt_horizons.to_csv( os.path.join(dir_res, "pred_horizons_alt.csv"),    index=False)
    rmst_alt_horizons.to_csv( os.path.join(dir_res, "rmst_horizons_alt.csv"),    index=False)

    print(f"    done replication {j}")

print(f"\nAll replications complete — results in results/sim_setting{SETTING}/{MODEL_NAME}/")
