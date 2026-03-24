# GAM-Based Delay Correction Model for Disease Surveillance

## Overview
This project develops a **Generalised Additive Model (GAM)-based framework** for correcting reporting delays in infectious disease surveillance data, with an application to **dengue incidence in Rio de Janeiro (2010–2013)**.

Real-time surveillance data is inherently incomplete due to reporting lags, leading to systematic underestimation of current incidence. This project addresses this by combining:

- **Nowcasting** (estimating unreported current cases)  
- **Forecasting** (predicting future incidence trends)  

The goal is to improve **real-time situational awareness and decision-making** under delayed data conditions.

---

## Key Features
- Delay correction using a **run-off triangle framework**
- Two GAM model specifications:
  - **Model Z**: Partial counts only (baseline)
  - **Model XZ**: Joint model using partial counts and truncated totals
- Captures **non-linear temporal and delay effects** using spline-based smooth terms
- Uses **Negative Binomial distribution** for overdispersed count data
- End-to-end pipeline:
  - Data preprocessing and masking (real-time simulation)
  - Model fitting (mgcv)
  - Predictive simulation (Monte Carlo)
  - Evaluation and visualisation

---

## Methodology

### Data Structure
- Weekly dengue case data organised as a **run-off triangle**
- Indexed by:
  - Time (epidemiological week)
  - Reporting delay (weeks)
- Missing values represent **unreported cases (nowcasting target)**

### Modelling Approach

Both models use GAMs to capture:
- Temporal trends  
- Reporting delay effects  
- Seasonality (week-of-year)

#### Model Z (Baseline)
- Models **partial counts only**
- Includes interaction between time and delay

#### Model XZ (Joint Model)
- Models both:
  - Partial counts
  - Truncated totals
- Improves explanatory power but introduces calibration trade-offs

---

## Evaluation Metrics

Models are evaluated across:

### Model Fit
- AIC  
- Deviance explained  

### Predictive Accuracy
- RMSE  
- MAE  

### Uncertainty
- 95% prediction interval coverage  
- Prediction interval width  

### Key Insight
- **Model Z** → more reliable uncertainty (conservative, wider intervals)  
- **Model XZ** → better fit and lower average error, but underestimates uncertainty  

---

## Results
- GAMs effectively capture **non-linear reporting delay structures**
- Joint modelling improves **average prediction error (MAE)**
- However, it reduces **uncertainty calibration**
- Highlights a key trade-off between:
  - Precision  
  - Reliability  

This trade-off is critical in **policy and operational settings**, where underestimating uncertainty can lead to poor decision-making.

---

## How to Run

### Requirements
- R (≥ 4.0)
- Required packages:
```r
data.table
mgcv
ggplot2
mvtnorm
dplyr
```

### Steps
1. Load dataset (.RData)
2. Run the main script:
```rsource("code_for_appendix.R")```
3. Outputs include:
- Model evaluation metrics
- Diagnostic plots
- Prediction interval comparisons
- Applications

---

Although developed for dengue surveillance, this framework is transferable to:

- Epidemiological nowcasting (e.g. COVID-19, influenza)
- Financial reporting delays
- ESG and climate data latency
- Any time-series with delayed observations

## Key Takeaways
GAMs provide a strong balance between:
- Flexibility
- Interpretability
- Computational efficiency

They offer a practical alternative to:
- Rigid GLMs
- Computationally intensive Bayesian models

---

## Future Work
- Incorporate spatial modelling (geographical variation)
- Extend to multi-disease applications
- Combine with Bayesian frameworks for improved uncertainty estimation
- Deploy as a real-time monitoring tool or dashboard

---
##Author

Linh Vu
MSc Environmental Intelligence, University of Exeter

---

##Acknowledgements

Supervised by Dr. Theo Economou
