# Reproduction materials: regression coefficients and cohort statistics from shared tree-ensemble models

Code, data and results for the paper *Regression coefficients and cohort statistics from shared tree-ensemble models without patient data: a simulation study and clinical applications*.

TreeIPF reads a shared tree-ensemble model file (rpart, XGBoost JSON, scikit-learn forests exported to JSON, ranger with `node.stats = TRUE`). From the trees and their leaf counts it rebuilds the joint distribution of the predictors by iterative proportional fitting, fits a linear or logistic working model to the ensemble's predictions inside the region covered by the splits, and takes the sample size from the root of the trees. A parametric bootstrap that retrains the ensemble gives a bias correction. The same functions measure how much of the cohort a file discloses.

## Contents

```
code/
  00_tree_ipf_lib.R          TreeIPF functions and the readers for the four file types
  sklearn_export.py          trains scikit-learn forests and saves their trees as JSON
  prepare_data.py            builds data/*.parquet from the raw IST and PhysioNet downloads
  01_run_simulations.R       Studies 1-5
  02_run_clinical_studies.R  Pima and PhysioNet audits (50 runs), PhysioNet XGBoost file, IST meta-analyses
  03_make_figures.R          figures
  07_run_bootstrap.R         Study 6, parametric bootstrap
  08_real_libraries.R        Study 7, XGBoost / scikit-learn / ranger files
  09_privacy_leakage.R       disclosure by tree depth and minimum leaf size
  10_summarise_results.R     performance measures with Monte Carlo standard errors
  11_make_tables.R           LaTeX tables from the result files
  12_checks.R, 12b_iterations_check.R   sensitivity checks (tail rule, grid error, residual variance, SMART updates)
  13_hyperparameters.R       Study 8, tuning and wrong hyperparameters
data/       analysis files, sources, licences and exclusion counts (see data/DATA_README.md)
results/    summary CSV files; replicates.zip holds every simulation replicate and clinical run
figures/    figure PDFs
```

## How to reproduce

R 4.4 with `rpart`, `data.table`, `Matrix`, `metafor`, `arrow`, `xgboost` (3.x), `ranger` (>= 0.16), `MASS`; Python 3 with `scikit-learn`, `pandas`, `pyarrow`. Run from the top folder of this repository (or set the environment variable `TREEIPF_ROOT`):

```
python code/prepare_data.py <folder with the raw IST and PhysioNet files>
Rscript code/01_run_simulations.R 1,2,3,4,5 8
Rscript code/02_run_clinical_studies.R 8
Rscript code/09_privacy_leakage.R 8
Rscript code/08_real_libraries.R 8
Rscript code/07_run_bootstrap.R 8 30 100
Rscript code/13_hyperparameters.R 8
Rscript code/12_checks.R 2
Rscript code/10_summarise_results.R
Rscript code/11_make_tables.R
Rscript code/03_make_figures.R
```

The last argument is the number of worker processes. Every replicate is saved as soon as it finishes and is skipped on a rerun; unzip `results/replicates.zip` into `results/reps/` to reuse ours. Seeds are fixed in the scripts. The analysis files in `data/` are already built, so `prepare_data.py` is needed only to rebuild them from the raw downloads.

## Data

- International Stroke Trial database (version 2), University of Edinburgh DataShare, https://doi.org/10.7488/ds/104, Open Data Commons Attribution License. Please cite Sandercock et al. (2011), Trials 12:101.
- PhysioNet/Computing in Cardiology Challenge 2012, https://physionet.org/content/challenge-2012/1.0.0/, Open Data Commons Attribution License. Please cite Silva et al. (2012) and Goldberger et al. (2000), Circulation 101:e215-e220.
- Pima Indians diabetes data from the R package `MASS` (Smith et al., 1988).

## Licence

Code: MIT (see LICENSE). The data keep their original licences.
