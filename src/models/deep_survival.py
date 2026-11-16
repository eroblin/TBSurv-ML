# =============================================================================
# src/models/deep_survival.py
# CoxCC and CoxTime models: building, hyperparameter search, and prediction.
# =============================================================================

import numpy as np
import pandas as pd
import optuna
import torch
import torchtuples as tt

from pycox.models import CoxTime, CoxCC
from pycox.models.cox_time import MLPVanillaCoxTime
from torchtuples.practical import MLPVanilla
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn_pandas import DataFrameMapper

# =============================================================================
# Activation and optimizer helpers
# =============================================================================

def get_activation(name):
    return {"relu": torch.nn.ReLU, "tanh": torch.nn.Tanh}[name]

def get_optimizer(name):
    return {
        "adam":         tt.optim.Adam,
        "adam_amsgrad": tt.optim.Adam,
        "rmsprop":      tt.optim.RMSprop,
        "sgdwr":        tt.optim.SGD,
    }[name]

# =============================================================================
# Model builder
# =============================================================================

def build_model(config, in_features, labtrans=""):
    out_features = labtrans.out_features if labtrans != "" else 1
    nb_neurons   = [config.neurons] * config.layers
    callbacks    = [tt.cb.EarlyStopping(load_best=False, checkpoint_model=False)]
    act          = get_activation(config.acti_func)
    optim        = get_optimizer(config.optim)(lr=config.lr, weight_decay=config.pen_l2)

    if config.name == "CoxCC":
        net   = MLPVanilla(in_features, nb_neurons, out_features,
                           batch_norm=True, dropout=config.dr,
                           activation=act, output_bias=False)
        model = CoxCC(net, optim)
    elif config.name == "CoxTime":
        net   = MLPVanillaCoxTime(in_features, nb_neurons,
                                  batch_norm=True, dropout=config.dr,
                                  activation=act)
        model = CoxTime(net, optim, labtrans=labtrans)

    if config.optim == "adam_amsgrad":
        optim.param_groups[0]["amsgrad"] = True
    elif config.optim == "sgdwr":
        optim.param_groups[0]["momentum"] = 0.9
        callbacks.append(tt.cb.LRCosineAnnealing())

    return model, callbacks

# =============================================================================
# Data mapper
# =============================================================================

def make_mapper(df):
    cont_vars = [c for c in df.columns if c not in ["yy", "status", "Trt"]]
    standardize = [([col], StandardScaler()) for col in cont_vars]
    leave       = [(col, None) for col in ["yy", "status", "Trt"]]
    return DataFrameMapper(standardize + leave, df_out=True)

# =============================================================================
# Hyperparameter search (Optuna objective)
# =============================================================================

def objective(trial, df, df_mapper, config):
    config.acti_func  = trial.suggest_categorical("activation",    ["relu", "tanh"])
    config.batch_size = trial.suggest_categorical("batch_size",    [8, 16, 32, 64, 128])
    config.layers     = trial.suggest_int("n_layers",              1, 4)
    config.lr         = trial.suggest_float("learning_rate",       1e-3, 1e-2)
    config.neurons    = trial.suggest_int("neurons",               4, 128)
    config.optim      = trial.suggest_categorical("optimizer",     ["adam", "adam_amsgrad", "rmsprop", "sgdwr"])
    config.pen_l2     = trial.suggest_float("l2",                  0.0, 0.1)
    config.dr         = trial.suggest_float("dropout",             0.0, 0.3)

    labtrans = CoxTime.label_transform() if config.name == "CoxTime" else ""

    kf        = KFold(n_splits=5, shuffle=True, random_state=397)
    val_loss  = pd.DataFrame()

    for t, (train_idx, val_idx) in enumerate(kf.split(df)):
        df_train = df_mapper.fit_transform(df.iloc[train_idx])
        df_val   = df_mapper.transform(df.iloc[val_idx]).astype("float32")

        x_train = np.array(df_train.drop(["yy", "status"], axis=1)).astype("float32")
        x_val   = np.array(df_val.drop(["yy", "status"],   axis=1)).astype("float32")
        y_train = (df_train["yy"].values, df_train["status"].values)
        y_val   = (df_val["yy"].values,   df_val["status"].values)

        if labtrans != "":
            y_train = labtrans.fit_transform(*y_train)
            y_val   = labtrans.transform(*y_val)

        val = tt.tuplefy(x_val, y_val).repeat(10).cat()

        model, callbacks = build_model(config, x_train.shape[1], labtrans)
        log = model.fit(x_train, y_train, int(config.batch_size),
                        epochs=500, callbacks=callbacks,
                        verbose=False, val_data=val)
        val_loss[t] = log.to_pandas()["val_loss"]

    return float(np.mean(val_loss).min())

# =============================================================================
# Fit final model with best hyperparameters
# =============================================================================

def fit_model(config, df_train, df_val, df_mapper):
    labtrans = CoxTime.label_transform() if config.name == "CoxTime" else ""

    df_train = df_mapper.fit_transform(df_train).astype("float32")
    df_val   = df_mapper.transform(df_val).astype("float32")

    x_train = np.array(df_train.drop(["yy", "status"], axis=1)).astype("float32")
    x_val   = np.array(df_val.drop(["yy", "status"],   axis=1)).astype("float32")
    y_train = (df_train["yy"].values, df_train["status"].values)
    y_val   = (df_val["yy"].values,   df_val["status"].values)

    if labtrans != "":
        y_train = labtrans.fit_transform(*y_train)
        y_val   = labtrans.transform(*y_val)

    val = tt.tuplefy(x_val, y_val).repeat(10).cat()

    model, callbacks = build_model(config, x_train.shape[1], labtrans)
    model.fit(x_train, y_train, int(config.batch_size),
              epochs=500, callbacks=callbacks,
              verbose=False, val_data=val)

    if config.name in ["CoxCC", "CoxTime"]:
        model.compute_baseline_hazards()

    return model, df_mapper

# =============================================================================
# Prediction
# =============================================================================

def predict_surv(model, df_test, df_mapper):
    df_test = df_mapper.transform(df_test).astype("float32")
    x_test  = np.array(df_test.drop(["yy", "status"], axis=1)).astype("float32")
    return model.predict_surv_df(x_test)