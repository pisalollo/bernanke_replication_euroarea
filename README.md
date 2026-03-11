# FAVAR Analysis of Monetary Policy Transmission in the Euro Area

## Project Overview
This repository contains a full-scale replication and extension of the **Factor-Augmented Vector Autoregression (FAVAR)** methodology proposed by **Bernanke, Boivin, and Eliasz (2005)**. 

The project applies this framework to the **Euro Area**, using a high-dimensional dataset to investigate how ECB monetary policy shocks propagate through a wide array of macroeconomic variables. Unlike standard VAR models, which are limited by the "curse of dimensionality," the FAVAR approach allows for a more comprehensive information set, mitigating the "price puzzle" and providing a richer description of the transmission mechanism.

## Key Features
- **High-Dimensional Dataset:** Processing and standardization of 100+ macroeconomic time series (Monthly, 2000-2025).
- **Factor Extraction & Identification:** - Dimensionality reduction via **Principal Component Analysis (PCA)**.
    - Implementation of **"Bernanke Cleaning"** to ensure the separation of latent factors from observed policy variables, dividing "Slow" and "Fast" variables.
    - Recursive identification (**Cholesky Decomposition**) and support for **Sign Restrictions**(to be implemented).
- **Monetary Policy Nuances:** Inclusion of the **Shadow Policy Rate** to account for the Zero Lower Bound (ZLB) period and unconventional monetary policy.
- **Robustness & Diagnostics:** - Automated lag selection via AIC/BIC.
    - Bootstrap simulations for IRF confidence intervals.
    - $R^2$ analysis and **Forecast Error Variance Decomposition (FEVD)**.
    - Historical Decomposition (to be implemented)

## Methodology

### 1. Data Preprocessing
The model classifies variables into **Fast** (e.g., exchange rates, financial indicators) and **Slow** moving (e.g., GDP, HICP) to correctly identify the policy shock, following the original Bernanke et al. (2005) identification scheme.

### 2. The Model
The economy is represented by the following transition equation:
$$Z_t = \Phi(L)Z_{t-1} + v_t$$
Where $Z_t$ contains both the extracted latent factors $F_t$ and the observed variables $Y_t$ (GDP, Inflation, and Policy Rate).

### 3. Impulse Response Functions (IRFs)
The code computes IRFs for both the variables included in the VAR and the "loaded" variables from the large dataset, allowing to observe the impact of a 25bps shock on indicators like:
- Real GDP & HICP
- Unemployment rate
- M1/M2 Monetary Aggregates
- Exchange Rates and Sentiment Indices (CCI/BCI)

## Repository Structure
- `main`: The primary script to run the analysis, from data loading to plot generation.
- `bernanke_cleaning.m` or 'bernanke_rotation.m': Function for extracting factors orthogonal to the policy rate.
- `CholeskyIdentification.m`: Core econometric routine for structural shock identification.
- `optimal_r.m`: Selects the optimal number of factors based on variance explained.

## How to Run
1. Ensure you have **MATLAB** installed with the *Econometrics Toolbox*.
2. Clone the repository.
3. Run `main`.

## Results
The model successfully identifies the transmission of monetary policy, showing a hump-shaped response of output and a gradual decline in price levels following a contractionary shock, consistent with economic theory and ECB's historical data.

---
**Author:** Lorenzo Pisa  
**Keywords:** Macroeconometrics, FAVAR, PCA, Monetary Policy, Euro Area, MATLAB.
