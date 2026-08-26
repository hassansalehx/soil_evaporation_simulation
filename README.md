# Soil evaporation simulation (SEC model)

Simulation of post-rainfall soil evaporation using the soil evaporation capacitance (SEC) model of [Or and Lehman (2019)](https://doi.org/10.1029/2018WR024050), applied to tropical-cyclone infiltration scenarios in the Arabian Peninsula.

This repository accompanies **Intensifying tropical cyclone in the Arabian Sea replenish depleting aquifers** ([Saleh et al., 2025](https://doi.org/10.1038/s43247-025-02493-w), *Communications Earth & Environment*).

**Author:** Hassan Saleh (Western Michigan University)

## Repository layout

```text
notebooks/   # soil_evap_sec.ipynb — main analysis
data/        # soil parameters + infiltration shapefiles (2011, 2018, 2020)
```

## Setup

From the repository root:

```bash
cd /path/to/this/repo
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -m ipykernel install --user --name soil-evap --display-name "Python (soil-evap)"
```

## Run

Open [`notebooks/soil_evap_sec.ipynb`](notebooks/soil_evap_sec.ipynb) with the **repository root as the working directory** and run cells sequentially.

The notebook loads soil parameters from `data/soil_param.xlsx` (from [Lehman et al., 2018](http://dx.doi.org/10.1029/2018GL078803)) and initial moisture shapefiles under `data/Infiltration_output/<year>/` for four soil types (Alluvium, Dunes, Hard Limestone, Karstified Limestone).

## Data

See [`data/README.md`](data/README.md) for file descriptions.

## Citation

Please cite the paper:

> Saleh, H., Sultan, M., Becker, R. et al. Intensifying tropical cyclone in the Arabian Sea replenish depleting aquifers. *Commun Earth Environ* **6**, 259 (2025). https://doi.org/10.1038/s43247-025-02493-w

Software metadata: [`CITATION.cff`](CITATION.cff).

## License

MIT — see [`LICENSE`](LICENSE).
