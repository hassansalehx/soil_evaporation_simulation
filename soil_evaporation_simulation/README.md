# Soil Evaporation Simulation

## Introduction
This notebook demonstrates the simulation of soil evaporation using the soil evaporation capacitance (SEC) model proposed by [Or and Lehman (2019)](https://doi.org/10.1029/2018WR024050). The goal is to simulate how soil moisture evaporates over time following tropical cyclone rainfall under different conditions using a capacitance approach. The model is initialized with moisture content (θ) and, based on soil parameters and the potential evaporation, the evaporation rate changes.

This notebook imports soil parameters based on the selected soil types that represent the study area taken from [Lehman et al. (2018)](http://dx.doi.org/10.1029/2018GL078803) and shapefiles containing the initial water content (in mm). The code processes the data according to the equations detailed in the supplementary materials of [Or and Lehman (2019)](https://doi.org/10.1029/2018WR024050) and visualizes the accumulated evaporation volume for each soil type.

This work is part of the study titled **"Intensifying Tropical Cyclone in the Arabian Sea Replenish Depleting Aquifers"** (under review in *Nature Communications Earth and Environment*).

The simulation proceeds by iteratively reducing the soil water content until the moisture reaches a residual value, applying different equations in:
- **Stage 1:** When water content is above the critical threshold.
- **Stage 2:** When water content falls below the critical threshold.

## Repository Structure

```
soil_evaporation_capacitance/
├── data/
│   ├── soil_param_related.xlsx
│   ├── Infiltration_output/
│   │   ├── 2011/
│   │   ├── 2018/
│   │   └── 2020/
├── notebooks/
│   └── soil_evap_sec.ipynb
├── README.md
├── requirements.txt
└── LICENSE



## Instructions
- **Data Organization:** Ensure all data files are placed in their respective folders as shown in the structure above.
- **Running the Notebook:** Open `soil_evap_sec.ipynb` in Jupyter and run the cells sequentially to load the data, perform the simulation, and view the plots.
- **Dependencies:** Install the required Python packages using:

pip install -r requirements.txt


## License
This project is licensed under the terms specified in the LICENSE file.

