# Data

`python code/prepare_data.py [raw-data folder]` builds `ist_dd6.parquet` and `physionet_ihd.parquet` from the raw downloads and writes the exclusion counts to `data_provenance.csv`. The rebuilt files are identical, value for value, to those used for the results.

## Sources and licences

**International Stroke Trial (IST)**, corrected database version 2, Edinburgh DataShare, doi:10.7488/ds/104, file `IST_corrected.csv` (19,435 patients). Licence: Open Data Commons Attribution (ODC-By), as in the licence file shipped with the dataset. Cite the dataset and Sandercock et al. (2011), *Trials* 12:101.

**PhysioNet/Computing in Cardiology Challenge 2012**, version 1.0.0, sets A, B and C (4,000 ICU stays each) with `Outcomes-a/b/c.txt`. Licence: Open Data Commons Attribution v1.0. Cite Goldberger et al. (2000), *Circulation* 101(23):e215-e220, and Silva et al. (2012), *Computing in Cardiology* 39:245-248.

## Outcomes

- IST, `y`: dead or dependent six months after randomisation (OCCODE 1 or 2) versus not recovered or recovered (3 or 4). OCCODE 0 or 9 (status unknown) is excluded.
- PhysioNet, `y`: in-hospital death (`In-hospital_death` = 1).

## IST variables

`id` is the row in `IST_corrected.csv` (from 0); `country` is CNTRYNUM (names in `ist_country_codes.csv`). Predictors, recorded at randomisation: `age` (years); `female` (1 = female); `rsbp` (systolic blood pressure, mmHg); `rconsc` (F fully alert, D drowsy, U unconscious); `ratrial` (atrial fibrillation, 1 = yes); `ct` (`none` = no CT before randomisation, `infarct` = CT showed a visible infarct, `no_infarct` = CT without a visible infarct); `rdef1` (face deficit) and `rdef2` (arm/hand deficit), 1 = yes, 0 = no or could not be assessed.

## PhysioNet variables

`id` is the RecordID. `age` (years) and `female` (1 = female; Gender 0 in the challenge coding) are admission descriptors. Over the first 48 hours of the ICU stay:

- `hr_max`: highest heart rate (beats/min).
- `sbp_min`: lowest systolic blood pressure, invasive or cuff (mmHg).
- `temp_max`: highest temperature (degrees C).
- `gcs_min`: lowest Glasgow Coma Score (3-15).
- `bun_max`: highest blood urea nitrogen (mg/dL).
- `creat_max`: highest serum creatinine (mg/dL).
- `wbc_max`: highest white blood cell count (10^3 cells/uL).
- `mechvent`: 1 if mechanical ventilation was charted at any time, else 0.

Implausible values are set to missing first (heart rate and blood pressure outside 1-300, temperature 25-45, GCS 3-15, BUN 0-300, creatinine 0.05-30, WBC 0.05-200).

## Exclusions (from `data_provenance.csv`)

IST: 19,435 randomised; 19,285 with known outcome; 18,304 complete on the eight predictors (all 981 lack atrial fibrillation, not recorded for the 984 pilot-phase patients); 15,580 in the 12 largest countries, which the analysis script selects.

PhysioNet: 12,000 stays, all with known outcome; 11,607 complete on the ten predictors, which the analysis script selects.

## Notes

- PhysioNet 2012 is a single source: MIMIC-II ICU stays from one hospital (sets A-C are random thirds of it).
- The IST groups are the countries of one randomised trial, not separate registries.
- Treatment allocation (aspirin, heparin) is not used. It was randomised, hence independent of the covariates, and the reference logistic regression and the tree method use exactly the same covariates, so the method comparison is like with like.
