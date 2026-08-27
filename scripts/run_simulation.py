#!/usr/bin/env python3
"""Run SEC evaporation simulation and save the manuscript-style figure."""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.sec_evaporation import (
    DEFAULT_YEARS,
    load_infiltration_gdfs,
    load_soil_params,
    plot_soil_evaporation,
    prepare_evaporation_data,
    repo_root,
    run_all_simulations,
)


def main() -> None:
    root = repo_root(ROOT)
    print(f"Repository root: {root}")

    soil_params = load_soil_params(root)
    gdfs = load_infiltration_gdfs(DEFAULT_YEARS, root)
    results = run_all_simulations(gdfs, soil_params, DEFAULT_YEARS)
    processed = prepare_evaporation_data(results)

    out = root / "outputs" / "figures" / "soil_evaporation_by_year.jpeg"
    plot_soil_evaporation(processed, save_path=out, show=False, verbose=True)
    print("Simulation complete.")


if __name__ == "__main__":
    main()
