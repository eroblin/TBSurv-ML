# TBSurv-ML

**R and Python code for evaluating prediction models that estimate individualized treatment benefit on time-to-event outcomes.**

This is the accompanying code of the paper:

> Elvire Roblin, Paul-Henry Cournède, Stefan Michiels  
> *Evaluation of machine-learning models to measure individualized treatment effects from randomized clinical trial data with time-to-event outcomes*  
> Cancer Informatics, 2026. DOI: 10.1177/11769351261491717  
> arXiv: https://arxiv.org/abs/XXXX.XXXXX

---

## Overview

The repository contains scripts to:

1. **Generate simulated data** for two data-generation processes:
   - **Setting 1** — Weibull survival times with nonlinear biomarker-by-treatment interactions (Biospear model)
   - **Setting 2** — Accelerated failure time (AFT) model with Friedman's random function generator

2. **Fit models** on the simulated data:
   - ALASSO (adaptive LASSO Cox proportional hazards, benchmark)
   - Interaction Forest (IF)
   - Causal Survival Forest (CSF)
   - CoxCC and CoxTime (feedforward neural networks, Python)

3. **Evaluate all models** and output the full set of performance measures:
   - Discrimination: C-statistic, C-for-Benefit
   - Calibration: E50-for-Benefit, ICI-for-Benefit
   - Prediction error: RMSE (PEHE)
   - Decision-making utility: AUTOC, Population Benefit (PB)

---

## Workflow

### Step 1 — Generate data

Open `scripts/01_generate_data.R`, set `SETTING <- 1` or `SETTING <- 2` at the top,
then source the file. Data are written under `data/sim_setting<k>/`.

### Step 2 — Fit a model (R)

Open `scripts/02_run_models.R`, set `SETTING` and `MODEL` at the top
(`"ALASSO"`, `"InteractionForest"`, or `"CSF"`), then source the file.
Predictions are written under `results/sim_setting<k>/<MODEL>/`.

### Step 2 — Fit a neural network (Python)

```bash
# edit model_name and setting at the top of the script, then:
python scripts/run_nn.py
```

### Step 3 — Evaluate

Open `scripts/03_evaluate.R`, set `SETTING` and `MODEL` to match step 2,
then source the file. All evaluation metrics are written to
`results/sim_setting<k>/<MODEL>/resAll_CI_match<MATCH_TYPE>.csv`.

---

## Simulation parameters (paper values)

### Setting 1

| Parameter | Value | Description |
|-----------|-------|-------------|
| n | 2 000 or 20 000 | total cohort size |
| p | 10 | number of active biomarkers |
| q | 20 | total number of biomarkers |
| α_τ | ln(0.75) ≈ −0.288 | treatment main effect |
| β_i1 | 0.288 | prognostic coefficient |
| β_i2 | −0.9 | interaction coefficient |
| f(X) | (2X − 1)² | nonlinear transformation |
| Censoring | 20 % or 50 % | independent |
| Replications | 100 | |

### Setting 2

| Parameter | Value | Description |
|-----------|-------|-------------|
| n | 2 000 or 20 000 | total cohort size |
| a_min / a_max | −0.2 / 0.3 | range for interaction coefficients a₂ₗ |
| n_tot | 5 | size of each random subset Rₗ |
| n_act | 5 (= l) | number of Friedman terms |
| l | 10 | number of terms in m₁(X) and m₂(X) |
| Censoring | 20 % | independent |
| Replications | 100 | |

---

## Dependencies

**R** (≥ 4.2): `survival`, `glmnet`, `ranger`, `grf`, `here`

**Python** (≥ 3.9): `pycox`, `torch`, `optuna`, `pandas`, `numpy`
