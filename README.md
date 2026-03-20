# FAVAR Analysis of Monetary Policy Transmission in the Euro Area

> Replication and extension of **Bernanke, Boivin & Eliasz (2005)** applied to the Euro Area using a high-dimensional macroeconomic dataset.

---

## Overview

This project implements a **Factor-Augmented Vector Autoregression (FAVAR)** model to investigate how ECB monetary policy shocks propagate through the Euro Area economy.

Standard VAR models are constrained by the *curse of dimensionality* and are prone to the *price puzzle*. The FAVAR framework addresses both issues by compressing a large dataset of 100+ macroeconomic indicators into a small set of latent factors, which are then included alongside explicit observable variables (GDP, HICP, Policy Rate) in the VAR system.

The project estimates five model variants and compares their impulse response functions:

| Model | Description |
|---|---|
| **VAR** | Baseline 3-variable VAR (GDP, HICP, Policy Rate) |
| **FAVAR** | Factors + GDP + HICP + Policy Rate |
| **FAVAR Baseline** | Factors + Policy Rate only |
| **FAVAR Shadow** | Factors + GDP + HICP + Shadow Rate |
| **FAVAR Baseline Shadow** | Factors + Shadow Rate only |

---

## Dataset

- **Source:** `data/EAdata2.xlsx`: https://sites.google.com/view/jingcynthiawu/shadow-rates
- **Coverage:** Monthly, January 2000 – November 2025 (~311 observations)
- **Dimensionality:** ~119 macroeconomic time series for the Euro Area
- **Explicit variables (Y block):** Real GDP (col. 1), HICP (col. 98), 3-Month Policy Rate (col. 75)
- **Wu-Xia Shadow Rate\*:** `data/shadowrate_spliced.mat` — (spliced after August 2022) used to handle the Zero Lower Bound (ZLB) period and unconventional monetary policy: https://sites.google.com/view/jingcynthiawu/shadow-rates

### Fast vs. Slow Variables

Variables are classified following Bernanke et al. (2005):

- **Fast-moving** (respond contemporaneously to policy): exchange rates (`ERUS`, `REER`), interest rates (`IRT`), equity indices (`SHIX`), sentiment (`BCI`, `CCI`, `SENTIX`), monetary aggregates (`M1_`, `M2_`, `CURR`), financial flows (`TASS`, `TLB`)
- **Slow-moving** (do not respond within the period): all remaining variables (GDP, HICP, unemployment, consumption, etc.)

This split is critical for the Bernanke rotation: slow variables are used to extract factors that are orthogonal to the policy instrument.

---

## Methodology

### 1. Factor Extraction — Bernanke Cleaning

Latent factors $\hat{F}_t$ are extracted via PCA and then *cleaned* to ensure they are orthogonal to the policy rate at impact. The number of factors $r$ is selected as the minimum needed to explain at least **60% of the total variance** (scree plot + cumulative variance criterion).

The transition equation of the FAVAR is:

$$Z_t = \Phi(L)\, Z_{t-1} + v_t, \qquad Z_t = [\hat{F}_t',\; Y_t']'$$

where $Y_t = [\text{GDP}_t,\; \text{HICP}_t,\; \text{Rate}_t]'$.

### 2. VAR Estimation and Lag Selection

- Lag order $p$ selected via **AIC** over $p \in \{1, \ldots, 12\}$
- Default: forced to $p = 7$ to capture monthly seasonality
- Residual autocorrelation is checked graphically after estimation

### 3. Identification — Cholesky Decomposition

Structural shocks are identified via **Cholesky decomposition** of the residual covariance matrix. The ordering places slow variables first and the policy rate last, so the rate shock has no contemporaneous effect on slow-moving variables.

### 4. IRF Computation and Bootstrap

- Horizon: **60 months**
- Confidence bands: **90% bootstrap intervals** (1000 replications, percentiles 5–95)
- IRFs are normalized to a **25 bps policy shock** via `normalize_irfToTarget`
- Factor IRFs are mapped back to observed variables via the **loading matrix** $\Lambda$ (`irf_loading`)

### 5. Diagnostics

- **R²** of the factor model for each variable in X (`calc_favar_r2`)
- **Forecast Error Variance Decomposition (FEVD)** at horizon 60 (`calc_favar_fevd`)
- **Robustness:** IRF comparison across $r = 1, \ldots, r_\text{opt}$ factors

---

## Repository Structure

```
├── main.m                      # Main script: data → estimation → IRFs → diagnostics
│
├── CholeskyIdentification.m# VAR estimation, Cholesky IRFs, bootstrap
├── optimal_r.m             # Scree plot + optimal number of factors (variance criterion)
├── get_indices.m           # Maps variable tickers to column indices in the dataset
├── calc_favar_r2.m         # R² of factor model for each observed variable
├── calc_favar_fevd.m       # Forecast Error Variance Decomposition
│
├── bernanke_cleaning/
│   ├── bernanke_cleaning.m     # Bernanke rotation with intercept (Fhat orthogonal to rate)
│   └── bernanke_rotation.m     # Bernanke rotation without intercept (alternative)
├── irf tools/
│   ├── irf_loading.m           # Projects factor IRFs onto observed variables via Lambda
│   ├── normalize_irfToTarget.m # Rescales IRFs to a target shock size (e.g. 25 bps)
│   ├── cumulative_irf.m        # Cumulates IRFs for selected variables (e.g. GDP, HICP)
│   ├── displayirf.m            # Plots IRF subpanels with bootstrap bands
│   └── compare_irf.m           # Overlays IRFs from multiple model variants
│
└── data/
    ├── EAdata2.xlsx            # Main dataset (100+ Euro Area monthly series)
    └── shadowrate_spliced.mat  # Shadow policy rate (ZLB adjustment)
```

> **Note:** `VarOLS`, `VarOLSbestP`, and `VarStr` are utility functions (not included above) required by `CholeskyIdentification` and `main.m`.

---

## How to Run

**Requirements:** MATLAB with the Econometrics Toolbox.


Key parameters at the top of `main.m`:

| Parameter | Default | Description |
|---|---|---|
| `confidence` | `[5 95]` | Bootstrap confidence level (90%) |
| `n_bootstrap` | `1000` | Number of bootstrap replications |
| `forcing` | `1` | Force a common lag length across all models |
| `forced_p` | `7` | Forced lag order (set to `0` for AIC selection) |
| `min_var_explained` | `0.60` | Minimum variance threshold for factor selection |

---

## Key Results

Following a 25 bps contractionary monetary policy shock:

- **Output (GDP):** hump-shaped contraction, peaking around 12–18 months
- **Prices (HICP):** gradual decline, no price puzzle
- **Unemployment:** rises with a lag
- **Monetary aggregates (M1, M2):** contract persistently
- **Exchange rate:** appreciates on impact
- **Consumer confidence (CCI):** immediate decline

Results are consistent with standard monetary transmission theory and ECB historical evidence. The FAVAR substantially outperforms the plain VAR in terms of R² and the absence of the price puzzle.

---

## Work in Progress

- [ ] **Historical Decomposition** — code is drafted and commented out at the bottom of `main.m`
- [ ] **Sign Restrictions** identification (alternative to Cholesky)

---

## Note

**Note on the Euro Area Shadow Rate**

This project uses the Wu & Xia (2016) shadow rate for the Euro Area
as a proxy for the monetary policy stance beyond the zero lower bound.
However, the EA shadow rate exhibits values reaching approximately -8%
during the negative interest rate policy (NIRP) period (2014–2022),
while the 3-month Euribor bottomed at around -0.5%.

This divergence raises concerns about the interpretability of the
shadow rate in the EA context. Unlike the US case — where the shadow
rate remains within a plausible range (-3% at trough) and cleanly
converges to the Fed funds rate upon lift-off — the EA shadow rate
appears to reflect model-implied extrapolation rather than an
economically meaningful measure of stance during the NIRP regime.

Following the authors' own recommendation, the spliced series
(shadow rate pre-ZLB exit, Euribor post-ZLB exit) is used in the
baseline specification. Results should be interpreted with caution. Alternative specifications using the
Euribor throughout are reported in the robustness checks.

**Current status: alternative measures under review**

---

**Author:** Lorenzo Pisa  
**Keywords:** Macroeconometrics · FAVAR · PCA · Monetary Policy · Euro Area · ECB · Shadow Rate · MATLAB