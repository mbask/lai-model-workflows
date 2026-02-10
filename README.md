# Reusable Bayesian Models and Workflows for LAI Estimation

This repository provides a reference implementation of a model- and workflow-based framework for Leaf Area Index (LAI) estimation in poplar systems. It accompanies the article submitted for publication and is intended as the first contribution to a shared repository of reusable, uncertainty-aware LAI estimation models.

The repository includes validation data, fitted Bayesian models, and fully reproducible R workflows designed to facilitate local recalibration and operational application across other forest contexts.

---

## Scope and Objectives

The primary objectives of this repository are to:

- Demonstrate a reproducible workflow for LAI estimation based on Bayesian hierarchical models;
- Provide validated model templates (AML and PAS) that can be reused and adapted to new sites;
- Promote a shift from static vegetation index–LAI parameterizations toward adaptable, data-driven modelling approaches;
- Support transparent uncertainty quantification in operational forest monitoring.

---

## Repository Contents

The repository is organized as follows:

  models/
    AML_PAS.rds        # Fitted AML and PAS models using NDVI and EVI Vegetation Indexes

  data/
    validation_LAI.csv  # Independent validation dataset

  scripts/
    validate_PAS.R   # Computation of RMSE, MAE, bias, and coverage on the validation dataset

  docs/
    AML_PAS_readme.md    # Model structure, prior specification, and validation dataset description

---

## Intended Use and Extension

This repository is provided as a reference implementation rather than a fixed solution. Users are encouraged to:

- Refit models using local calibration data;
- Adapt prior distributions to site-specific knowledge;
- Extend the workflow to additional forest types or vegetation indices;
- Contribute new validated models following the same structure.

---

## Citation

If you use this repository, please cite both the associated article and the Zenodo archive:
(to be completed)

> Author(s). Year. *Reusable Models and Workflows for LAI Estimation*. Zenodo. DOI: XX.XXXX/zenodo.XXXXXXX

---

## License

This repository is released under a MIT license to facilitate reuse and extension. See the LICENSE file for details.



