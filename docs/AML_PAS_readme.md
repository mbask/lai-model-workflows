# AML-PAS Bayesian models documentation

## Validation Dataset

### Description of the Validation Data

The validation dataset used in this study consists of independent LAI measurements not employed during model fitting. LAI observations were derived from ground-based daily measurements  collected in hybrid poplar plantations (Chianucci, F., Bergante, S., Chiarabaglio, P. M., & Bascietto, M. (2025). Monitoring hybrid poplar plantations using continuous canopy photography: Influence of clone and water status. Annals of Silvicultural Research, 50(2), 52–57. https://doi.org/10.12899/asr-2679) and aggregated at weekly temporal resolution.

All validation data refer to irrigated stands of the *I-214* poplar clone with a stand age of five years, ensuring independence from the training dataset while maintaining ecological comparability. 

Validation data are provided as a serialized objects file `data/validation_dataset.rds` holding a list of two datasets:

    > summary(vld_l)
         Length Class      Mode
    EVI  7      data.table list
    NDVI 7      data.table list

Datasets are stored as a `data.frame/data.table` to ensure efficient manipulation and reproducibility. Sample size (N) of the `data.table`s is 90 observations. Each row corresponds to a  daily plot-level observation, combining field-measured LAI, remotely sensed vegetation indices.

The datasets include:

- A plot identifier (`plotID`, only for reference to
  [@chianucci_monitoring_2025], character class)

- Leaf Area Index observations (`LAI_r`, numeric class);

- Corresponding vegetation indices values (`index_value`, scaled in
  0--10,000 range, numeric class);

- Type of vegetation index (`index`, coded as factor: \"NDVI\",
  \"EVI\");

- Stand age `i.e.` years since establishment (`Age`, numeric class);

- Week of year (`week`, numeric class);

- Phenological month (`phenoMonth`, coded as factor: \"06\", \"07\",
  \"09\");

The validation dataset can be loaded in `R` environment as:

      vld_l <- readRDS("data/validation_dataset.rds")

The `scripts/validate_PAS.R` script provides a fully reproducible workflow to
validate the VI-based Bayesian LAI models reported in the manuscript and
to reproduce Figure 4. It loads the independent validation dataset
(`validation_dataset.rds`) and the fitted `brms` model objects and
training data (`model.rds`). Because models were trained on standardized
predictors and a standardized LAI response, the script (i) rescales
validation predictors (Age and vegetation index values) using the
training-set centering and scaling parameters, and (ii) back-transforms
posterior predictions to LAI units. For each vegetation index (NDVI,
EVI), population-level posterior expected predictions are generated via
`brms::posterior_epred(..., re_formula = NA)`, thereby excluding
group-level effects to assess transferability to new plots. Predicted
summaries (posterior mean and 95% credible interval) are then joined to
the validation observations and aggregated to weekly means. Finally, the
script computes global and phenology-stratified validation metrics
(RMSE, MAE, bias, empirical 95% coverage, and interval width) and
produces the multi-panel agreement/coverage plots composing Figure 4.


## Model Specifications

The fitted model objects (AML and PAS) are `brmsfit` instances. It
contains: (i) the full model specification (formula, family, link
functions), (ii) prior definitions, (iii) data and preprocessing
metadata, (iv) MCMC configuration and sampling output (posterior draws),
and (v) convergence and sampling diagnostics (R-hat, ESS, divergent
transitions). All posterior summaries (fixed effects, smooth terms,
random effects) and posterior predictive quantities reported here were
derived directly from this object to ensure full reproducibility.

All 4 models (AML, PAS models based on EVI and NDVI VI datasets) are
provided as serialized `brmsfit` objects, enabling direct reuse for
inspection or comparison, in file `data/model.rds`.

The Age--Monotonic Logistic (AML) models describe the ontogenetic
trajectory of LAI as a function of stand age using a sigmoidal
formulation. Separate AML models were fitted for NDVI- and EVI-based
predictors.

Model structure:

- Formula: `LAI ~ index_value * Age + (1 | phenoMnth)`

- Gaussian likelihood for LAI;

- Nonlinear logistic mean function;

- Fixed effects only;

- Weakly informative priors on asymptote, slope, and inflection point.

Phenology-Adaptive Spline (PAS) models extend AML formulations by
allowing the vegetation index--LAI relationship to vary across
phenological months.

Key features:

- Formula:
  `LAI ~ index_value + s(Age, k = 5) + (1 + index_value | phenoMnth)`

- Smooth age effect with an estimated breakpoint representing canopy
  closure;

- Random intercepts and slopes by phenological month;

- Bayesian hierarchical structure enabling partial pooling;

- Explicit uncertainty propagation for derived quantities (e.g., canopy
  closure age).

### Priors

Prior distributions were chosen to be weakly to moderately informative
and are summarized in Table [1].

:::
| Model |              Parameter              |     Prior      |       Rationale       |
|-------|-------------------------------------|----------------|-----------------------|
| AML   | Intercept                           | Normal(0, 1.5) | Centered LAI          |
|       | index_value`                        | Normal(0, 0.7) | Plausible effect size |
|       | `Age`                               | Normal(0, 0.5) | Moderate age effect   |
|       | `index_value` $\times$ `Age`        | Normal(0, 0.3) | Strong shrinkage      |
|       | `phenoMnth` SD                      | Exponential(2) | Few levels            |
|       | Residual SD ($\sigma$)              | Exponential(1) | Weakly informative    |
| PAS   | Intercept                           | Normal(0, 1.5) | Centered LAI          |
|       | `index_value`                       | Normal(0, 0.7) | Plausible slope       |
|       | Smooth `Age` SD                     | Exponential(2) | Penalized smooth      |
|       | `phenoMnth` SD (intercept)          | Exponential(2) | Few levels            |
|       | `phenoMnth` SD (`index_value`slope) | Exponential(3) | Strong shrinkage      |
|       | Random-effect correlation           | LKJ(3)         | Weak correlations     |
|       | Residual SD ($\sigma$)              | Exponential(1) | Weakly informative    |


  : Prior distributions used for the Bayesian AML and PAS models. All
  predictors and the LAI response were standardized prior to model
  fitting.
:::

### EVI-based models

Details of AML and PAS models trained on EVI data are found in
Table [2].

:::
|                | `fit_evi_l$AML` | `fit_evi_l$PAS` |
|----------------|-----------------|-----------------|
| Family         | gaussian        | gaussian        |
| Link           | identity        | identity        |
| Algorithm      | sampling        | sampling        |
| Chains         | 4               | 4               |
| Iter           | 5000            | 5000            |
| Warmup         | 1500            | 1500            |
| N              | 51              | 51              |
| Rhat_max_fixed | 1.001184        | 1.000896        |
| ESS_min_fixed  | 4397.545        | 5091.830        |
| Divergences    | 0               | 0               |
| Treedepth_hits | 0               | 0               |

  : Model specification and MCMC diagnostics for the EVI-based AML and
  PAS Bayesian models fitted with brms, including sampling settings,
  sample size, and convergence indicators (R-hat, ESS, and NUTS
  diagnostics).
:::

### NDVI-based models

Details of AML and PAS models trained on EVI data are found in
Table [3].

:::

|                | `fit_ndvi_l$AML` | `fit_ndvi_l$PAS` |
|----------------|------------------|------------------|
| Family         | gaussian         | gaussian         |
| Link           | identity         | identity         |
| Algorithm      | sampling         | sampling         |
| Chains         | 4                | 4                |
| Iter           | 5000             | 5000             |
| Warmup         | 1500             | 1500             |
| N              | 51               | 51               |
| Rhat_max_fixed | 1.0005491.000560 |                  |
| ESS_min_fixed  | 4383.170         | 6189.077         |
| Divergences    | 0                | 0                |
| Treedepth_hits | 0                | 0                |


  : Model specification and MCMC diagnostics for the NDVI-based AML and
  PAS Bayesian models fitted with brms, including sampling settings,
  sample size, and convergence indicators (R-hat, ESS, and NUTS
  diagnostics).
:::

## R Scripts and Workflow

### Model Application Script

An annotated R script (`scripts/validate_PAS.R`) is provided to demonstrate how the
PAS model can be applied to new data. The script performs the following
steps:

1.  Data import and formatting;

2.  Predictor standardization using training-data statistics;

3.  Posterior prediction of LAI;

4.  Extraction of posterior means and credible intervals;

5.  Calculation of validation metrics (RMSE, MAE, bias, coverage).

The script is designed to be modular and easily adaptable to other
sites, forest types, or vegetation indices.

### Software and Dependencies

All analyses were conducted in `R` using the following main packages:

- `brms` for Bayesian modelling (required);

- `data.table` for data handling; (required to run the script but
  fallback to base `data.frame` data structure is possible)

- `ggplot2`, and `patchwork` for plotting manuscript figure (optional)

Package versions and session information are included at the end of the
script to ensure full reproducibility.
