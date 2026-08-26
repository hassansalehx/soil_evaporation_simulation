# Data

## `soil_param.xlsx`

Soil hydraulic parameters for the four study soil types, used by the SEC model. Values follow [Lehman et al. (2018)](http://dx.doi.org/10.1029/2018GL078803); the notebook reads sheet `code`.

## `Infiltration_output/<year>/`

Initial post-infiltration soil moisture (mm) by soil type, exported as ESRI shapefiles for three tropical-cyclone years:

| Year | Event context (see paper) |
|------|---------------------------|
| 2011 | — |
| 2018 | — |
| 2020 | — |

Each year folder contains shapefiles for:

- Alluvium
- Dunes
- Hard Limestone
- Karstified Limestone

These inputs drive the evaporation simulations in `notebooks/soil_evap_sec.ipynb`.
