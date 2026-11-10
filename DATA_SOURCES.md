# Data sources and attribution

The fitted-cost acquisition experiments use the following public datasets. Their original raw archives are not duplicated in this repository; `acquisition/scripts/download_uci_data.py` retrieves them from UCI. The frozen processed windows, derived interval-cost boxes, and experiment outputs are included under the repository's CC BY 4.0 data license.

1. **Occupancy Detection** — Luis Candanedo (2016), UCI Machine Learning Repository, DOI: [10.24432/C5X01N](https://doi.org/10.24432/C5X01N). UCI lists the dataset under CC BY 4.0. Download: https://archive.ics.uci.edu/static/public/357/occupancy%2Bdetection.zip
2. **EEG Eye State** — Oliver Roesler (2013), UCI Machine Learning Repository, DOI: [10.24432/C57G7J](https://doi.org/10.24432/C57G7J). UCI lists the dataset under CC BY 4.0. Download: https://archive.ics.uci.edu/static/public/264/eeg%2Beye%2Bstate.zip
3. **SECOM** — Michael McCann and Adrian Johnston (2008), UCI Machine Learning Repository, DOI: [10.24432/C54305](https://doi.org/10.24432/C54305). UCI lists the dataset under CC BY 4.0. Download: https://archive.ics.uci.edu/static/public/179/secom.zip

The selected windows and preprocessing rules are documented in `acquisition/docs/ACQUISITION_REPRODUCTION.md`. The source pipeline preserves the original feature ordering, stable time ordering, imputation, standardization, intercept, and linear time-trend construction.

The SKAB n=768 stress case is cited separately and not redistributed because its upstream repository/Kaggle distribution carries separate copyleft terms. See `docs/SKAB_REPRODUCTION.md`.
